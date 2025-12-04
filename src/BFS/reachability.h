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

  const Graph &G;
  size_t round;
  hashbag<NodeId> bag;
  sequence<NodeId> frontier;
  sequence<uint8_t> visited;
  sequence<uint8_t> in_frontier;
  constexpr static bool use_local_queue = true;

#ifdef STAT
  size_t max_frontier_size;
  size_t min_frontier_size;
  size_t total_frontier_size;
  size_t total_visited_vertices_in_queue;
  size_t total_visited_edges_in_queue;
#endif
 public:
  size_t LOCAL_QUEUE_SIZE = 32;
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
#ifdef STAT
    parlay::sequence<size_t> visited_vertices(frontier_size);
    parlay::sequence<size_t> visited_edges(frontier_size);
#endif
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
#ifdef STAT
                visited_vertices[i]++;
                visited_edges[i] += deg;
#endif
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
#ifdef STAT
    total_visited_vertices_in_queue += parlay::reduce(visited_vertices);
    total_visited_edges_in_queue += parlay::reduce(visited_edges);
    total_frontier_size += frontier_size;
    double avg_visited_vertices =
        1.0 * parlay::reduce(visited_vertices) / frontier_size;
    size_t max_visited_vertices =
        parlay::reduce(visited_vertices, maxm<size_t>());
    size_t min_visited_vertices =
        parlay::reduce(visited_vertices, minm<size_t>());
    double avg_visited_edges =
        1.0 * parlay::reduce(visited_edges) / frontier_size;
    size_t max_visited_edges = parlay::reduce(visited_edges, maxm<size_t>());
    size_t min_visited_edges = parlay::reduce(visited_edges, minm<size_t>());
    cout << "round " << round
         << ": avg_visited_vertices: " << avg_visited_vertices
         << ", max_visited_vertices: " << max_visited_vertices
         << ", min_visited_vertices: " << min_visited_vertices << endl;
    cout << "\tavg_visited_edges: " << avg_visited_edges
         << ", max_visited_edges: " << max_visited_edges
         << ", min_visited_edges: " << min_visited_edges << endl;
#endif
  }

  sequence<uint8_t> reachability(NodeId s) {
#ifdef STAT
    parlay::internal::timer tm;
    max_frontier_size = 0;
    min_frontier_size = G.n;
    total_frontier_size = 0;
    total_visited_vertices_in_queue = 0;
    total_visited_edges_in_queue = 0;
#endif

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
#ifdef STAT
      max_frontier_size = max(max_frontier_size, frontier_size);
      min_frontier_size = min(min_frontier_size, frontier_size);
#endif
      printf("Round %zu: size: %zu, local: %d, ", round, frontier_size,
             use_local_queue);
      relax(frontier_size);
      // t.next("sparse");
      round++;
    }
#ifdef STAT
    tm.stop();
    size_t total_visited_vertices = parlay::count(visited, true);
    double total_time = tm.total_time();
    std::cout << "LOCAL_QUEUE_SIZE: " << LOCAL_QUEUE_SIZE << '\n';
    std::cout << "max_frontier_size: " << max_frontier_size << '\n';
    std::cout << "min_frontier_size: " << min_frontier_size << '\n';
    std::cout << "total_frontier_size: " << total_frontier_size << '\n';
    std::cout << "total_visited_vertices_in_queue: "
              << total_visited_vertices_in_queue << '\n';
    std::cout << "total_visited_edges_in_queue: "
              << total_visited_edges_in_queue << '\n';
    std::cout << "total_visited_vertices: " << total_visited_vertices << '\n';
    std::cout << "total_rounds: " << round << '\n';
    std::cout << "time: " << total_time << '\n';
    ofstream ofs("reachability_stats.tsv", ios::app);
    ofs << fixed << setprecision(6);
    ofs << LOCAL_QUEUE_SIZE << '\t' << max_frontier_size << '\t'
        << min_frontier_size << '\t' << total_frontier_size << '\t'
        << total_visited_vertices_in_queue << '\t'
        << total_visited_edges_in_queue << '\t' << total_visited_vertices
        << '\t' << round << '\t' << total_time << '\n';
    ofs.close();
#endif
    return visited;
  }
};
