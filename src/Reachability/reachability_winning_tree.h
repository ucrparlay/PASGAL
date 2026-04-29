#pragma once
#include <algorithm>
#include <atomic>
#include <climits>
#include <vector>

#include "graph.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "winning_tree.h"

using namespace std;
using namespace parlay;

// Reachability with the WinningTree (64-ary hierarchical bitmap) as the
// frontier collector, following the BFS/bfs_delta_winning_tree.h pattern.
//
// Phase concurrency: WinningTree supports concurrent inserts OR iteration,
// not both.  We use two trees (curr / next) and swap each round so iteration
// of curr happens while inserts go to next.
template <class Graph>
class ReachabilityWinningTree {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr size_t CACHELINE_SIZE = 64;
  static constexpr size_t EDGE_PER_CACHELINE =
      CACHELINE_SIZE / sizeof(typename Graph::Edge);

  const Graph &G;
  WinningTree curr_bag_;
  WinningTree next_bag_;
  sequence<atomic<bool>> visited;
  size_t beta_;
  size_t max_queue_size_;
  constexpr static bool use_local_queue = true;

 public:
  ReachabilityWinningTree() = delete;
  ReachabilityWinningTree(const Graph &_G, size_t beta = 2048,
                          size_t max_queue_size = 1000)
      : G(_G),
        curr_bag_(G.n),
        next_bag_(G.n),
        beta_(beta),
        max_queue_size_(max_queue_size) {
    visited = sequence<atomic<bool>>(G.n);
  }

  size_t beta() const { return beta_; }
  size_t max_queue_size() const { return max_queue_size_; }

  // No in_frontier dedup needed: WinningTree::insert is idempotent (it just
  // sets a bit), and the visited CAS upstream guarantees at most one
  // discoverer per (round, vertex) anyway.
  void add_to_frontier(NodeId v) { next_bag_.insert(v); }

  void visit_neighbors_parallel(NodeId u) {
    parallel_for(
        G.offsets[u], G.offsets[u + 1],
        [&](size_t i) {
          NodeId v = G.edges[i].v;
          if (!visited[v].exchange(true)) {
            add_to_frontier(v);
          }
        },
        beta_);
  }

  void visit_neighbors_sequential(NodeId u, NodeId *local_queue, size_t &rear) {
    for (EdgeId i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      NodeId v = G.edges[i].v;
      if (!visited[v].exchange(true)) {
        if (rear < max_queue_size_) {
          local_queue[rear++] = v;
        } else {
          add_to_frontier(v);
        }
      }
    }
  }

  void relax() {
    [[maybe_unused]] static const int num_threads = parlay::num_workers();
    size_t frontier_size = std::max((size_t)1, curr_bag_.size());
    const size_t queue_size = std::min(
        max_queue_size_,
        std::max((size_t)1, num_threads * beta_ / frontier_size));

    curr_bag_.iterate_all_and_clear([&](size_t fv) {
      NodeId f = static_cast<NodeId>(fv);
      if constexpr (use_local_queue) {
        thread_local std::vector<NodeId> local_queue_storage;
        if (local_queue_storage.size() < max_queue_size_) {
          local_queue_storage.resize(max_queue_size_);
        }
        NodeId *local_queue = local_queue_storage.data();
        size_t front = 0, rear = 0;
        size_t vertices_visited = 0;
        size_t edges_processed = 0;
        const size_t max_edges = queue_size * EDGE_PER_CACHELINE;
        local_queue[rear++] = f;
        while (front < rear && vertices_visited < max_queue_size_ &&
               edges_processed < max_edges) {
          NodeId u = local_queue[front++];
          vertices_visited++;
          size_t deg = G.offsets[u + 1] - G.offsets[u];
          edges_processed += deg;
          if (deg < beta_) {
            visit_neighbors_sequential(u, local_queue, rear);
          } else {
            visit_neighbors_parallel(u);
          }
        }
        while (front < rear) {
          add_to_frontier(local_queue[front++]);
        }
      } else {
        visit_neighbors_parallel(f);
      }
    });
  }

  sequence<bool> reachability(NodeId s) {
    parallel_for(0, G.n, [&](size_t i) {
      visited[i].store(false, std::memory_order_relaxed);
    });

    visited[s].store(true, std::memory_order_relaxed);
    curr_bag_.insert(s);

    while (true) {
      if (curr_bag_.empty()) break;
      relax();
      std::swap(curr_bag_, next_bag_);
    }
    return parlay::tabulate(G.n, [&](size_t i) -> bool {
      return visited[i].load(std::memory_order_relaxed);
    });
  }
};
