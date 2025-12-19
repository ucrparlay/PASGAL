#pragma once
#include <atomic>
#include <climits>
#include <functional>

#include "graph.h"
#include "hashbag.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"

using namespace std;
using namespace parlay;

template <class Graph>
class BFS {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr NodeId DIST_MAX = numeric_limits<NodeId>::max();
  static constexpr size_t LOCAL_QUEUE_SIZE = 128;
  static constexpr size_t BLOCK_SIZE = 1024;
  static constexpr size_t NUM_SAMPLES = 1024;
  static constexpr size_t SPARSE_TH = 20;
  static constexpr size_t GROWTH_FACTOR = 10;

  const Graph &G;
  int step;
  int bag_id;
  const size_t num_bags;
  size_t round;
  sequence<hashbag<NodeId>> bags;
  sequence<NodeId> frontier;
  sequence<NodeId> dist;
  sequence<int> steps;
  sequence<std::atomic<bool>> in_frontier;
  bool sparse;
  constexpr static bool use_local_queue = true;

  template <class T, class Less = std::less<T>>
  static bool write_min_std(T *addr, T value, Less less = Less()) {
    std::atomic_ref<T> ref(*addr);
    T current = ref.load(std::memory_order_acquire);
    while (less(value, current)) {
      if (ref.compare_exchange_weak(current, value, std::memory_order_acq_rel,
                                    std::memory_order_relaxed)) {
        return true;
      }
    }
    return false;
  }

 public:
  BFS() = delete;
  BFS(const Graph &_G, int _step)
      : G(_G), step(_step), num_bags((LOCAL_QUEUE_SIZE + step - 1) / step + 2) {
    bags = sequence<hashbag<NodeId>>(num_bags, hashbag<NodeId>(G.n));
    frontier = sequence<NodeId>::uninitialized(G.n);
    dist = sequence<NodeId>::uninitialized(G.n);
    steps = sequence<int>(G.n, INT_MAX);
    in_frontier = sequence<std::atomic<bool>>::uninitialized(G.n);
  }

  void add_to_frontier(NodeId v) {
    bool expected = false;
    in_frontier[v].compare_exchange_strong(
        expected, true, std::memory_order_acq_rel, std::memory_order_acquire);
    if (sparse) {
      // Add v to the next bag even if it belongs to the current one
      int id = max<int>(bag_id + 1, dist[v] / step);
      assert(id - bag_id <= (int)num_bags);
      if (write_min_std(&steps[v], id)) {
        bags[id % num_bags].insert(v);
      }
    }
  }

  size_t estimate_size() {
    static uint32_t seed = 951;
    size_t hits = 0;
    for (size_t i = 0; i < NUM_SAMPLES; i++) {
      NodeId u = hash32(seed) % G.n;
      if (dist[u] == round) {
        hits++;
      }
      seed++;
    }
    return hits * G.n / NUM_SAMPLES;
  }

  void visit_neighbors_parallel(NodeId u) {
    parallel_for(
        G.offsets[u], G.offsets[u + 1],
        [&](size_t i) {
          NodeId v = G.edges[i].v;
          if (write_min_std(&dist[v], dist[u] + 1)) {
            add_to_frontier(v);
          }
        },
        BLOCK_SIZE);
  }

  void visit_neighbors_sequential(NodeId u, NodeId *local_queue, size_t &rear) {
    for (EdgeId i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      NodeId v = G.edges[i].v;
      if (write_min_std(&dist[v], dist[u] + 1)) {
        if (rear < LOCAL_QUEUE_SIZE) {
          local_queue[rear++] = v;
        } else {
          add_to_frontier(v);
        }
      }
    }
  }

  void dense2sparse() {
    for (size_t i = 0; i < num_bags; i++) {
      bags[i].clear();
    }
    parallel_for(0, G.n, [&](size_t i) {
      if (in_frontier[i].load(std::memory_order_acquire)) {
        bags[dist[i] % num_bags].insert(i);
      }
    });
  }

  void sparse_relax(size_t frontier_size) {
    parallel_for(0, frontier_size, [&](size_t i) {
      NodeId f = frontier[i];
      in_frontier[f].store(false, std::memory_order_release);
      if (bag_id == steps[f]) {
        steps[f] = INT_MAX;
        if (use_local_queue) {
          NodeId local_queue[LOCAL_QUEUE_SIZE];
          size_t front = 0, rear = 0;
          local_queue[rear++] = f;
          while (front < rear) {
            NodeId u = local_queue[front++];
            size_t deg = G.offsets[u + 1] - G.offsets[u];
            if (deg < BLOCK_SIZE) {
              visit_neighbors_sequential(u, local_queue, rear);
            } else {
              visit_neighbors_parallel(u);
            }
          }
        } else {
          visit_neighbors_parallel(f);
        }
      }
    });
  }

  void dense_relax() {
    parallel_for(0, G.n, [&](NodeId u) {
      if (dist[u] > round + 1) {
        const auto neighbors = G.in_neighors(u);
        for (size_t j = 0; j < neighbors.size(); j++) {
          NodeId v = neighbors[j].v;
          if (dist[v] != DIST_MAX && dist[u] > dist[v] + 1) {
            dist[u] = dist[v] + 1;
            in_frontier[u].store(true, std::memory_order_release);
            if (dist[v] == round) {
              break;
            }
          }
        }
      } else if (dist[u] <= round) {
        if (in_frontier[u].load(std::memory_order_acquire)) {
          in_frontier[u].store(false, std::memory_order_release);
        }
      }
    });
  }

  bool if_sparse(size_t frontier_size) {
    return frontier_size * SPARSE_TH < G.n;
  }

  sequence<NodeId> bfs(NodeId s) {
    parallel_for(0, G.n, [&](size_t i) {
      in_frontier[i].store(false, std::memory_order_release);
      dist[i] = DIST_MAX;
    });

    sparse = true;
    dist[s] = 0;
    steps[s] = 0;
    in_frontier[s].store(true, std::memory_order_release);

    round = 0;
    bag_id = 0;
    bags[0].insert(s);

    int last_update = 0;
    while (true) {
      internal::timer t;
      size_t frontier_size =
          bags[bag_id % num_bags].pack_into(make_slice(frontier));
      if (frontier_size) {
        last_update = round;
      } else {
        if (round - last_update > num_bags) {
          break;
        }
      }
      printf("Round %zu: bag_id: %d, size: %zu, local: %d, ", round, bag_id,
             frontier_size, use_local_queue);
      sparse_relax(frontier_size);
      sparse = true;
      t.next("sparse");
      round++;
      bag_id++;
    }
    printf("final round: %zu\n", round);

#if 1
    for (size_t i = 0; i < num_bags; i++) {
      assert(bags[i].pack_into(make_slice(frontier)) == 0);
    }
#endif
    return dist;
  }
};
