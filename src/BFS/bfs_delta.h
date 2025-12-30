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
  static constexpr size_t BLOCK_SIZE = 1024;
  static constexpr size_t MAX_QUEUE_SIZE = 128;
  static constexpr size_t CACHELINE_SIZE = 64;
  static constexpr size_t EDGE_PER_CACHELINE =
      CACHELINE_SIZE / sizeof(typename Graph::Edge);

  const Graph &graph_;
  NodeId threshold_;
  NodeId delta_;
  hashbag<NodeId> curr_bag_;
  hashbag<NodeId> next_bag_;
  sequence<NodeId> frontier_;
  sequence<NodeId> dist_;
  sequence<std::atomic<bool>> in_frontier_;
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
  BFS(const Graph &graph, NodeId delta)
      : graph_(graph), delta_(delta), curr_bag_(graph.n), next_bag_(graph.n) {
    frontier_ = sequence<NodeId>::uninitialized(graph.n);
    dist_ = sequence<NodeId>::uninitialized(graph.n);
    in_frontier_ = sequence<std::atomic<bool>>::uninitialized(graph.n);
  }

  void add_to_frontier(NodeId v) {
    bool expected = false;
    if (in_frontier_[v].compare_exchange_strong(expected, true,
                                                std::memory_order_acq_rel,
                                                std::memory_order_acquire)) {
      if (dist_[v] <= threshold_) {
        curr_bag_.insert(v);
      } else {
        next_bag_.insert(v);
      }
    }
  }

  void visit_neighbors_parallel(NodeId u) {
    parallel_for(
        graph_.offsets[u], graph_.offsets[u + 1],
        [&](size_t i) {
          NodeId v = graph_.edges[i].v;
          if (write_min_std(&dist_[v], dist_[u] + 1)) {
            add_to_frontier(v);
          }
        },
        BLOCK_SIZE);
  }

  void visit_neighbors_sequential(NodeId u, NodeId *local_queue, size_t &rear) {
    for (EdgeId i = graph_.offsets[u]; i < graph_.offsets[u + 1]; i++) {
      NodeId v = graph_.edges[i].v;
      if (write_min_std(&dist_[v], dist_[u] + 1)) {
        if (rear < MAX_QUEUE_SIZE) {
          local_queue[rear++] = v;
        } else {
          add_to_frontier(v);
        }
      }
    }
  }

  void sparse_relax(size_t frontier_size) {
    parallel_for(0, frontier_size, [&](size_t i) {
      NodeId f = frontier_[i];
      in_frontier_[f].store(false, std::memory_order_release);
      if constexpr (use_local_queue) {
        NodeId local_queue[MAX_QUEUE_SIZE];
        size_t front = 0, rear = 0;
        size_t vertices_visited = 0;
        size_t edges_processed = 0;
        const size_t max_edges = MAX_QUEUE_SIZE * EDGE_PER_CACHELINE;
        local_queue[rear++] = f;
        while (front < rear && vertices_visited < MAX_QUEUE_SIZE &&
               edges_processed < max_edges) {
          NodeId u = local_queue[front++];
          vertices_visited++;
          size_t deg = graph_.offsets[u + 1] - graph_.offsets[u];
          edges_processed += deg;
          if (deg < BLOCK_SIZE) {
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

  sequence<NodeId> bfs(NodeId s) {
    parallel_for(0, graph_.n, [&](size_t i) {
      in_frontier_[i].store(false, std::memory_order_release);
      dist_[i] = DIST_MAX;
    });
    dist_[s] = 0;
    in_frontier_[s].store(true, std::memory_order_release);
    threshold_ = delta_;

    [[maybe_unused]] int round = 0;
    curr_bag_.insert(s);
    while (true) {
      internal::timer t;
      size_t frontier_size = curr_bag_.pack_into(make_slice(frontier_));
      if (frontier_size == 0) {
        break;
      }
      while (frontier_size) {
        sparse_relax(frontier_size);
        frontier_size = curr_bag_.pack_into(make_slice(frontier_));
      }
      round++;
      std::swap(curr_bag_, next_bag_);
      threshold_ += delta_;
    }
    return dist_;
  }
};
