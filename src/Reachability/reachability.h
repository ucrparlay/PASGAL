#pragma once
#include <alloca.h>

#include <atomic>
#include <climits>

#include "graph.h"
#include "hashbag.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"

using namespace std;
using namespace parlay;

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
  //                      threshold for visiting a single vertex's
  //                      neighbors.
  //   max_queue_size  -- vertex cap on each task's local queue.  Each
  //                      task drains up to this many vertices from the
  //                      bag before flushing remaining discoveries
  //                      back to the bag.
  //
  // The local queue can be very large (up to ~1M vertices on
  // high-diameter graphs), so it is heap-allocated per worker thread
  // (`thread_local`) rather than placed on the stack.  Each worker pays
  // the allocation once and reuses the buffer across tasks.
  //
  // Defaults below are from the 2026-04-29 sweep across 21 graphs on
  // sloth (2x Xeon Skylake, 96 HT): (beta=2048, q=2000) wins both
  // geomean (1.13x of per-graph oracle) and worst-case (1.70x).  Adaptive
  // threshold modes (graph_deg / frontier_deg, edge / dual) were
  // explored and dropped -- vertex with a constant cap dominated.  See
  // results/SWEEP_SUMMARY.md.
  // ------------------------------------------------------------------

  const Graph &G;
  hashbag<NodeId> bag;
  sequence<NodeId> frontier;
  sequence<atomic<bool>> visited;
  size_t beta_;
  size_t max_queue_size_;
  constexpr static bool use_local_queue = true;

 public:
  Reachability() = delete;
  Reachability(const Graph &_G, size_t beta = 2048,
               size_t max_queue_size = 2000)
      : G(_G), bag(G.n), beta_(beta), max_queue_size_(max_queue_size) {
    frontier = sequence<NodeId>::uninitialized(G.n);
    visited = sequence<atomic<bool>>(G.n);
  }

  size_t beta() const { return beta_; }
  size_t max_queue_size() const { return max_queue_size_; }

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

  void visit_neighbors_sequential(NodeId u, NodeId *local_queue, size_t &rear) {
    // Prefetch distance for visited[]: G.edges[] is scanned sequentially so
    // hardware prefetch keeps it warm; visited[] is fully random so we issue
    // explicit prefetches a fixed number of edges ahead.
    constexpr EdgeId EDGE_PF = 16;
    const EdgeId end = G.offsets[u + 1];
    for (EdgeId i = G.offsets[u]; i < end; i++) {
      if (i + EDGE_PF < end)
        __builtin_prefetch(&visited[G.edges[i + EDGE_PF].v], 0, 0);
      NodeId v = G.edges[i].v;
      if (!visited[v].exchange(true, std::memory_order_relaxed)) {
        if (rear < max_queue_size_) {
          // Prefetch v's offset so it is hot when v is later popped.
          __builtin_prefetch(&G.offsets[v], 0, 1);
          local_queue[rear++] = v;
        } else {
          add_to_frontier(v);
        }
      }
    }
  }

  void relax(size_t frontier_size) {
    parallel_for(
        0, frontier_size,
        [&](size_t i) {
          NodeId f = frontier[i];
          if constexpr (use_local_queue) {
            // Stack-allocated per task via alloca.  thread_local would
            // alias between sibling tasks (the same worker can run
            // other outer iterations while this one is suspended inside
            // visit_neighbors_parallel's nested parallel_for, clobbering
            // the shared buffer).  alloca lives on the worker's call
            // stack and is freed automatically when this lambda
            // invocation returns; standard C++ has no portable way to
            // size a stack array at runtime.  At
            // max_queue_size=2000 * sizeof(NodeId)=4 bytes, that's 8KB
            // per active task -- trivial for a Parlay worker stack.
            NodeId *local_queue = static_cast<NodeId *>(
                alloca(max_queue_size_ * sizeof(NodeId)));
            size_t front = 0, rear = 0;
            size_t vertices_visited = 0;
            local_queue[rear++] = f;
            // Burst: drain the local queue until we've processed
            // max_queue_size_ vertices or the queue is empty.  Anything
            // still in the queue at burst-end gets flushed to the bag.
            // Prefetch distance for the edge-list start: relies on
            // G.offsets[w] being warm (seeded by the push-time prefetch
            // in visit_neighbors_sequential above).
            constexpr size_t VERTEX_PF = 8;
            while (front < rear && vertices_visited < max_queue_size_) {
              if (front + VERTEX_PF < rear) {
                NodeId w = local_queue[front + VERTEX_PF];
                __builtin_prefetch(&G.edges[G.offsets[w]], 0, 0);
              }
              NodeId u = local_queue[front++];
              vertices_visited++;
              size_t deg = G.offsets[u + 1] - G.offsets[u];
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
