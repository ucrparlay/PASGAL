#pragma once
#include <climits>

#include "graph.h"
#include "hashbag.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "utils.h"

using namespace std;
using namespace parlay;

template <class Graph>
class BFS {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr NodeId DIST_MAX = numeric_limits<NodeId>::max();
  static constexpr size_t LOCAL_QUEUE_SIZE = 1;
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
  sequence<uint8_t> in_frontier;
  bool sparse;
  constexpr static bool use_local_queue = false;

 public:
  BFS() = delete;
  BFS(const Graph &_G, int _step)
      : G(_G), step(_step), num_bags((LOCAL_QUEUE_SIZE + step - 1) / step + 2) {
    bags = sequence<hashbag<NodeId>>(num_bags, hashbag<NodeId>(G.n));
    frontier = sequence<NodeId>::uninitialized(G.n);
    dist = sequence<NodeId>::uninitialized(G.n);
    steps = sequence<int>(G.n, INT_MAX);
    in_frontier = sequence<uint8_t>::uninitialized(G.n);
  }

  void add_to_frontier(NodeId v) {
    if (in_frontier[v] == false) {
      in_frontier[v] = true;
    }
    if (sparse) {
      // Add v to the next bag even if it belongs to the current one
      int id = max<int>(bag_id + 1, dist[v] / step);
      assert(id - bag_id <= (int)num_bags);
      if (write_min(&steps[v], id)) {
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
          if (write_min(&dist[v], dist[u] + 1)) {
            add_to_frontier(v);
          }
        },
        BLOCK_SIZE);
  }

  void visit_neighbors_sequential(NodeId u, NodeId *local_queue, size_t &rear) {
    for (EdgeId i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      NodeId v = G.edges[i].v;
      if (write_min(&dist[v], dist[u] + 1)) {
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
      if (in_frontier[i]) {
        bags[dist[i] % num_bags].insert(i);
      }
    });
  }

  void sparse_relax(size_t frontier_size) {
    parallel_for(0, frontier_size, [&](size_t i) {
      NodeId f = frontier[i];
      in_frontier[f] = false;
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
            in_frontier[u] = true;
            if (dist[v] == round) {
              break;
            }
          }
        }
      } else if (dist[u] <= round) {
        if (in_frontier[u]) {
          in_frontier[u] = false;
        }
      }
    });
  }

  bool if_sparse(size_t frontier_size) {
    return frontier_size * SPARSE_TH < G.n;
  }

  sequence<NodeId> bfs(NodeId s) {
    parallel_for(0, G.n, [&](size_t i) {
      in_frontier[i] = false;
      dist[i] = DIST_MAX;
    });

    sparse = true;
    dist[s] = 0;
    steps[s] = 0;
    in_frontier[s] = true;

    round = 0;
    bag_id = 0;
    bags[0].insert(s);

    // int last_update = 0;
    while (true) {
      size_t approx_size = estimate_size();
      if (if_sparse(approx_size)) {
        if (!sparse) {
          dense2sparse();
        }
        size_t frontier_size =
            bags[bag_id % num_bags].pack_into(make_slice(frontier));
        if (!frontier_size) {
          break;
        }
        // printf("Round %zu: size: %zu, local: %d, ", round, frontier_size,
        // use_local_queue);
        sparse_relax(frontier_size);
        sparse = true;
        // t.next("sparse");
      } else {
        // printf("Round %zu: ", round);
        dense_relax();
        sparse = false;
        // t.next("dense");
      }
      round++;
      bag_id++;
    }
    return dist;
  }
};
