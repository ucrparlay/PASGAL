#pragma once
#include <climits>

#include "chunkbag.h"
#include "graph.h"
#include "hashbag.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "utils.h"

using namespace std;
using namespace parlay;

template <class Graph, template <class> class BagT = hashbag>
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
  const int LOG2N;
  const size_t num_bags;
  size_t round;
  sequence<BagT<NodeId>> bags;
  sequence<NodeId> frontier;
  sequence<NodeId> dist;
  sequence<int> bag_id;
  sequence<atomic<bool>> in_frontier;
  bool sparse;
  bool use_local_queue;

  // dist is written concurrently (write_min in the sparse phase, owner
  // stores in the dense phase) while other tasks read it, so all
  // accesses inside parallel regions go through atomic_ref.
  NodeId load_dist(NodeId v) {
    return std::atomic_ref<NodeId>(dist[v]).load(std::memory_order_relaxed);
  }
  void store_dist(NodeId v, NodeId d) {
    std::atomic_ref<NodeId>(dist[v]).store(d, std::memory_order_relaxed);
  }

  // Bag sizing lives here, not in the bag: this frontier re-inserts a
  // vertex every time write_min shrinks its bag_id (and the flag test
  // above it is a plain load/store, so concurrent discoverers can both
  // insert), so unlike the delta-stepping BFS a bag can exceed n live
  // elements.  Keep the historical 2n bound, plus one partial chunk per
  // worker for chunk bags.
  static size_t bag_capacity(size_t n) {
    if constexpr (requires { BagT<NodeId>::CHUNK_SIZE; }) {
      return 2 * n + parlay::num_workers() * BagT<NodeId>::CHUNK_SIZE;
    } else {
      return n;  // hashbag sizes its own table from the element count
    }
  }

 public:
  BFS() = delete;
  BFS(const Graph &_G)
      : G(_G), LOG2N(log2_up(G.n)), num_bags(log2_up(LOCAL_QUEUE_SIZE) + 2) {
    bags = sequence<BagT<NodeId>>(num_bags, BagT<NodeId>(bag_capacity(G.n)));
    frontier = sequence<NodeId>::uninitialized(G.n);
    dist = sequence<NodeId>::uninitialized(G.n);
    bag_id = sequence<int>::uninitialized(G.n);
    in_frontier = sequence<atomic<bool>>(G.n);
  }

  // `first` is true when the caller knows v was just discovered
  // (write_min replaced DIST_MAX), which makes the in_frontier test
  // unnecessary: v cannot be in the frontier yet.
  void add_to_frontier(NodeId v, bool first = false) {
    if (sparse) {
      NodeId dist_v = load_dist(v);
      int id = dist_v == 0 ? 0 : log2_up(dist_v);
      if (first || in_frontier[v].load(std::memory_order_relaxed) == false) {
        in_frontier[v].store(true, std::memory_order_relaxed);
        write_min(&bag_id[v], id);
        bags[id % num_bags].insert(v);
      } else {
        if (write_min(&bag_id[v], id)) {
          bags[id % num_bags].insert(v);
        }
      }
    } else {
      in_frontier[v].store(true, std::memory_order_relaxed);
    }
  }

  size_t estimate_size([[maybe_unused]] size_t id) {
    static uint32_t seed = 951;
    size_t hits = 0;
    for (size_t i = 0; i < NUM_SAMPLES; i++) {
      NodeId u = hash32(seed) % G.n;
      if (load_dist(u) == round) {
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
          NodeId old_dist;
          if (write_min_old(&dist[v], load_dist(u) + 1, old_dist)) {
            add_to_frontier(v, old_dist == DIST_MAX);
          }
        },
        BLOCK_SIZE);
  }

  void visit_neighbors_sequential(NodeId u, NodeId *local_queue, size_t &rear) {
    for (EdgeId i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      NodeId v = G.edges[i].v;
      NodeId old_dist;
      if (write_min_old(&dist[v], load_dist(u) + 1, old_dist)) {
        if (rear < LOCAL_QUEUE_SIZE) {
          local_queue[rear++] = v;
        } else {
          add_to_frontier(v, old_dist == DIST_MAX);
        }
      }
    }
  }

  void dense2sparse() {
    for (size_t i = 0; i < num_bags; i++) {
      bags[i].clear();
    }
    parallel_for(0, G.n, [&](size_t i) {
      if (in_frontier[i].load(std::memory_order_relaxed)) {
        size_t id = log2_up(load_dist(i));
        bags[id % num_bags].insert(i);
      }
    });
  }

  void sparse_relax(size_t id, size_t frontier_size) {
    parallel_for(0, frontier_size, [&](size_t i) {
      NodeId f = frontier[i];
      in_frontier[f].store(false, std::memory_order_relaxed);
      if (id == 0 || id == log2_up(load_dist(f))) {
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

  void dense_relax([[maybe_unused]] size_t id) {
    parallel_for(0, G.n, [&](NodeId u) {
      // Only this task writes dist[u] in the dense phase, so du tracks
      // it exactly; neighbors' distances are racy reads by design.
      NodeId du = load_dist(u);
      if (du > round + 1) {
        const auto neighbors = G.in_neighors(u);
        for (size_t j = 0; j < neighbors.size(); j++) {
          NodeId v = neighbors[j].v;
          NodeId dv = load_dist(v);
          if (dv != DIST_MAX && du > dv + 1) {
            du = dv + 1;
            store_dist(u, du);
            in_frontier[u].store(true, std::memory_order_relaxed);
            if (dv == round) {
              break;
            }
          }
        }
      } else if (du <= round) {
        // Conditional store: unconditionally writing false would dirty
        // Theta(n) cache lines per dense round for settled vertices.
        if (in_frontier[u].load(std::memory_order_relaxed)) {
          in_frontier[u].store(false, std::memory_order_relaxed);
        }
      }
    });
  }

  bool if_sparse(size_t frontier_size) {
    return frontier_size * SPARSE_TH < G.n;
  }

  void prepare(NodeId s) {
    parallel_for(0, G.n, [&](size_t i) {
      in_frontier[i] = false;
      dist[i] = DIST_MAX;
      bag_id[i] = LOG2N;
    });

    sparse = true;
    dist[s] = 0;
    add_to_frontier(s);

    round = 0;
  }

  void run_prepared() {
    size_t prev_size = 0;
    bool dense = false;
    for (int i = 0; i <= LOG2N; i++) {
      if (i != 0) {
        round = max(round, (size_t)((1 << (i - 1)) + 1));
      }
      while (true) {
        size_t approx_size = estimate_size(i);
        // internal::timer t;
        // printf("prev_size: %zu, approx_size: %zu\n", prev_size, approx_size);
        if (if_sparse(approx_size)) {
          if (dense) {
            dense2sparse();
          }
          size_t frontier_size =
              bags[i % num_bags].pack_into(make_slice(frontier));
          if (frontier_size <= prev_size * GROWTH_FACTOR) {
            use_local_queue = true;
          } else {
            use_local_queue = false;
          }
          prev_size = frontier_size;
          if (!frontier_size) {
            break;
          }
          // printf("Round %zu: size: %zu, local: %d, ", round, frontier_size,
          // use_local_queue);
          sparse_relax(i, frontier_size);
          dense = false;
          round++;
          // t.next("sparse");
        } else {
          // printf("Round %zu: ", round);
          dense_relax(i);
          dense = true;
          round++;
          // t.next("dense");
          prev_size = approx_size;
        }
      }
    }

    for (size_t i = 0; i < num_bags; i++) {
      assert(bags[i].pack_into(make_slice(frontier)) == 0);
    }
  }

  const sequence<NodeId> &distances() const { return dist; }

  sequence<NodeId> bfs(NodeId s) {
    prepare(s);
    run_prepared();
    return dist;
  }
};
