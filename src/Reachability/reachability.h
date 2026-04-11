#pragma once
#include <climits>

#include "graph.h"
#include "hashbag.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "utils.h"

using namespace std;
using namespace parlay;

template <class Graph>
class Reachability {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr NodeId DIST_MAX = numeric_limits<NodeId>::max();
  static constexpr size_t BLOCK_SIZE = 1024;
  static constexpr size_t NUM_SAMPLES = 1024;
  static constexpr size_t SPARSE_TH = 20;
  static constexpr size_t LOCAL_QUEUE_SIZE = 5000;

  const Graph &G;
  size_t round;
  hashbag<NodeId> bag;
  sequence<NodeId> frontier;
  sequence<uint8_t> visited;
  sequence<uint8_t> in_frontier;
  constexpr static bool use_local_queue = true;

 public:
  Reachability() = delete;
  Reachability(const Graph &_G) : G(_G), bag(G.n) {
    frontier = sequence<NodeId>::uninitialized(G.n);
    visited = sequence<uint8_t>::uninitialized(G.n);
    in_frontier = sequence<uint8_t>::uninitialized(G.n);
  }

  void add_to_frontier(NodeId v) {
    if (compare_and_swap<uint8_t>(&in_frontier[v], false, true)) {
      bag.insert(v);
    }
  }

  void visit_neighbors_parallel(NodeId u) {
    parallel_for(G.offsets[u], G.offsets[u + 1], [&](size_t i) {
      NodeId v = G.edges[i].v;
      if (compare_and_swap<uint8_t>(&visited[v], false, true)) {
        add_to_frontier(v);
      }
    });
  }

  void visit_neighbors_sequential(NodeId u, NodeId *local_queue, size_t &rear) {
    for (EdgeId i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      NodeId v = G.edges[i].v;
      if (compare_and_swap<uint8_t>(&visited[v], false, true)) {
        if (rear < LOCAL_QUEUE_SIZE) {
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
          in_frontier[f] = false;
          if (use_local_queue) {
            NodeId local_queue[LOCAL_QUEUE_SIZE];
            size_t front = 0, rear = 0, sum = 0;
            local_queue[rear++] = f;
            while (front < rear && sum < LOCAL_QUEUE_SIZE) {
              NodeId u = local_queue[front++];
              size_t deg = G.offsets[u + 1] - G.offsets[u];
              if (deg < LOCAL_QUEUE_SIZE) {
                visit_neighbors_sequential(u, local_queue, rear);
                sum += deg;
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

  sequence<uint8_t> reachability(NodeId s) {
    parallel_for(0, G.n, [&](size_t i) {
      in_frontier[i] = false;
      visited[i] = false;
    });

    round = 0;
    visited[s] = true;
    add_to_frontier(s);

    while (true) {
      // internal::timer t;
      size_t frontier_size = bag.pack_into(make_slice(frontier));
      if (!frontier_size) {
        break;
      }
      relax(frontier_size);
      // t.next("sparse");
      round++;
    }
    return visited;
  }
};
