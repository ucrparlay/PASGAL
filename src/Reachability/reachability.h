#pragma once
#include <algorithm>
#include <atomic>
#include <climits>
#include <vector>

#include "graph.h"
#include "hashbag.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"

using namespace std;
using namespace parlay;

// Threshold for ending each per-task local-queue burst:
//   DUAL         = stop when EITHER vertices_visited >= max_queue_size OR
//                  edges_processed >= max_edges (current behaviour).
//   VERTEX       = stop when vertices_visited >= max_queue_size.  No edge
//                  bound; no adaptive scaling.
//   EDGE         = stop when edges_processed >= max_edges, where max_edges
//                  uses adaptive scaling
//                  (queue_size = min(max_queue_size, num_threads*beta/frontier_size)).
//                  No vertex bound.
//   GRAPH_DEG    = vertex bound, but cap = clamp(c_const / (G.m/G.n), 10,
//                  max_queue_size).  Computed once at construction.
//   FRONTIER_DEG = vertex bound, cap = clamp(c_const / frontier_avg_deg,
//                  10, max_queue_size).  Recomputed each relax() call
//                  using the average out-degree of the current frontier.
enum class ThresholdMode {
  DUAL = 0,
  VERTEX = 1,
  EDGE = 2,
  GRAPH_DEG = 3,
  FRONTIER_DEG = 4,
};

template <class Graph>
class Reachability {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  // ------------------------------------------------------------------
  // Tunable parameters.  Sweep via the CLI to find good values for a
  // given graph and machine.
  //
  //   beta            -- edges per burst.  Used as parallel_for grain
  //                      and as the deg-based sequential vs parallel
  //                      threshold.
  //   max_queue_size  -- vertex cap on each task's local queue.
  //                      Decoupled from beta because Reachability has
  //                      no level-correctness constraint (BFS does).
  //
  // The local queue can be very large (up to ~1M vertices on
  // high-diameter graphs), so it is heap-allocated per worker thread
  // (`thread_local`) rather than placed on the stack.  Each worker pays
  // the allocation once and reuses the buffer across tasks.
  // ------------------------------------------------------------------
  static constexpr size_t CACHELINE_SIZE = 64;
  static constexpr size_t EDGE_PER_CACHELINE =
      CACHELINE_SIZE / sizeof(typename Graph::Edge);

  static constexpr size_t Q_MIN = 1;

  const Graph &G;
  hashbag<NodeId> bag;
  sequence<NodeId> frontier;
  sequence<atomic<bool>> visited;
  size_t beta_;
  size_t max_queue_size_;
  size_t c_const_;          // numerator for GRAPH_DEG / FRONTIER_DEG
  size_t graph_deg_q_;      // precomputed cap for GRAPH_DEG
  ThresholdMode mode_;
  constexpr static bool use_local_queue = true;

 public:
  Reachability() = delete;
  // Defaults chosen from the 2026-04-28 sweep across 12 graphs, 5 sources
  // each: VERTEX-only threshold with q=2000 minimizes geomean runtime
  // and gives the widest "near-optimal" plateau (q in [1000, 20000]).
  Reachability(const Graph &_G, size_t beta = 2048,
               size_t max_queue_size = 2000,
               ThresholdMode mode = ThresholdMode::VERTEX,
               size_t c_const = 0)
      : G(_G), bag(G.n), beta_(beta), max_queue_size_(max_queue_size),
        c_const_(c_const), mode_(mode) {
    frontier = sequence<NodeId>::uninitialized(G.n);
    visited = sequence<atomic<bool>>(G.n);
    // Precompute the per-graph effective q for GRAPH_DEG mode.
    double avg_deg = G.n ? double(G.m) / double(G.n) : 1.0;
    if (avg_deg < 1.0) avg_deg = 1.0;
    size_t q = c_const_ ? size_t(double(c_const_) / avg_deg) : max_queue_size_;
    if (q < Q_MIN) q = Q_MIN;
    if (q > max_queue_size_) q = max_queue_size_;
    graph_deg_q_ = q;
  }

  size_t beta() const { return beta_; }
  size_t max_queue_size() const { return max_queue_size_; }
  size_t c_const() const { return c_const_; }
  ThresholdMode mode() const { return mode_; }
  const char *mode_str() const {
    switch (mode_) {
      case ThresholdMode::DUAL:         return "dual";
      case ThresholdMode::VERTEX:       return "vertex";
      case ThresholdMode::EDGE:         return "edge";
      case ThresholdMode::GRAPH_DEG:    return "graph_deg";
      case ThresholdMode::FRONTIER_DEG: return "frontier_deg";
    }
    return "?";
  }

  // Each caller of add_to_frontier(v) has just won
  // visited[v].exchange(true), so v is unique across all concurrent
  // workers in this round.  No dedup needed -- a plain bag.insert is
  // sufficient.
  void add_to_frontier(NodeId v) { bag.insert(v); }

  void visit_neighbors_parallel(NodeId u) {
    parallel_for(
        G.offsets[u], G.offsets[u + 1],
        [&](size_t i) {
          NodeId v = G.edges[i].v;
          // memory_order_relaxed is sufficient: visited[v] is just a
          // dedup flag, no happens-before relation needs to be
          // established with any other vertex's state.
          if (!visited[v].exchange(true, std::memory_order_relaxed)) {
            add_to_frontier(v);
          }
        },
        beta_);
  }

  void visit_neighbors_sequential(NodeId u, NodeId *local_queue, size_t &rear,
                                  size_t cap) {
    for (EdgeId i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      NodeId v = G.edges[i].v;
      if (!visited[v].exchange(true, std::memory_order_relaxed)) {
        if (rear < cap) {
          local_queue[rear++] = v;
        } else {
          add_to_frontier(v);
        }
      }
    }
  }

  // Compute the effective vertex cap for this relax() based on mode.
  size_t effective_q(size_t frontier_size) {
    switch (mode_) {
      case ThresholdMode::GRAPH_DEG:
        return graph_deg_q_;
      case ThresholdMode::FRONTIER_DEG: {
        // Average out-degree of the current frontier.
        size_t total_edges = parlay::reduce(parlay::delayed_tabulate(
            frontier_size, [&](size_t i) -> size_t {
              NodeId u = frontier[i];
              return G.offsets[u + 1] - G.offsets[u];
            }));
        double avg = frontier_size ? double(total_edges) / double(frontier_size)
                                   : 1.0;
        if (avg < 1.0) avg = 1.0;
        size_t q = c_const_ ? size_t(double(c_const_) / avg) : max_queue_size_;
        if (q < Q_MIN) q = Q_MIN;
        if (q > max_queue_size_) q = max_queue_size_;
        return q;
      }
      default:
        return max_queue_size_;
    }
  }

  void relax(size_t frontier_size) {
    [[maybe_unused]] static const int num_threads = parlay::num_workers();
    const size_t cap = effective_q(frontier_size);
    const size_t queue_size = std::min(
        cap,
        std::max((size_t)1, num_threads * beta_ / frontier_size));

    parallel_for(
        0, frontier_size,
        [&](size_t i) {
          NodeId f = frontier[i];
          if constexpr (use_local_queue) {
            // Heap-allocated, reused per worker.  Sized once when the
            // worker first sees a queue this big.
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
            // Hot loop: branch on mode_ once via a lambda so the compiler
            // can specialize.  mode_ is constant for the whole burst, so
            // the branch predictor learns it instantly anyway.
            auto keep_going = [&]() {
              if (front >= rear) return false;
              switch (mode_) {
                case ThresholdMode::DUAL:
                  return vertices_visited < cap &&
                         edges_processed < max_edges;
                case ThresholdMode::VERTEX:
                case ThresholdMode::GRAPH_DEG:
                case ThresholdMode::FRONTIER_DEG:
                  return vertices_visited < cap;
                case ThresholdMode::EDGE:
                  return edges_processed < max_edges;
              }
              return false;
            };
            while (keep_going()) {
              NodeId u = local_queue[front++];
              vertices_visited++;
              size_t deg = G.offsets[u + 1] - G.offsets[u];
              edges_processed += deg;
              if (deg < beta_) {
                visit_neighbors_sequential(u, local_queue, rear, cap);
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
        },
        1);
  }

  sequence<bool> reachability(NodeId s) {
    parallel_for(0, G.n, [&](size_t i) {
      visited[i].store(false, std::memory_order_relaxed);
    });

    visited[s].store(true, std::memory_order_relaxed);
    add_to_frontier(s);

    while (true) {
      size_t frontier_size = bag.pack_into(make_slice(frontier));
      if (!frontier_size) break;
      relax(frontier_size);
    }
    return parlay::tabulate(G.n, [&](size_t i) -> bool {
      return visited[i].load(std::memory_order_relaxed);
    });
  }
};
