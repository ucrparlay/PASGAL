#pragma once
#include <atomic>
#include <climits>
#include <functional>

#include "graph.h"
#include "hashbag.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/utilities.h"

using namespace std;
using namespace parlay;

template <class Graph>
class BFS {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr NodeId DIST_MAX = numeric_limits<NodeId>::max();
  static constexpr size_t BETA = 2048;
  static constexpr size_t CACHELINE_SIZE = 64;
  static constexpr size_t EDGE_PER_CACHELINE =
      CACHELINE_SIZE / sizeof(typename Graph::Edge);  // 16 for 32-bit edge
  static constexpr size_t MAX_QUEUE_SIZE = BETA / EDGE_PER_CACHELINE;  // 128

  const Graph &graph_;
  NodeId threshold_;
  NodeId delta_;
  hashbag<NodeId> curr_bag_;
  hashbag<NodeId> next_bag_;
  sequence<NodeId> frontier_;
  sequence<NodeId> dist_;
  sequence<std::atomic<bool>> in_curr_frontier_;
  sequence<std::atomic<bool>> in_next_frontier_;
  constexpr static bool use_local_queue = true;

  NodeId load_dist(NodeId v,
                   std::memory_order order = std::memory_order_acquire) {
    return std::atomic_ref<NodeId>(dist_[v]).load(order);
  }

  template <class T, class Less = std::less<T>>
  static bool write_min_std(T *addr, T value, Less less = Less()) {
    std::atomic_ref<T> ref(*addr);
    T current = ref.load(std::memory_order_acquire);
    while (less(value, current)) {
      if (ref.compare_exchange_weak(current, value, std::memory_order_acq_rel,
                                    std::memory_order_acquire)) {
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
    in_curr_frontier_ = sequence<std::atomic<bool>>::uninitialized(graph.n);
    in_next_frontier_ = sequence<std::atomic<bool>>::uninitialized(graph.n);
  }

  void add_to_frontier(NodeId v, NodeId dist_v) {
    bool expected = false;
    if (dist_v < threshold_) {
      if (in_curr_frontier_[v].compare_exchange_strong(
              expected, true, std::memory_order_acq_rel,
              std::memory_order_relaxed)) {
        curr_bag_.insert(v);
      }
    } else {
      assert(dist_v == threshold_);
      if (in_next_frontier_[v].compare_exchange_strong(
              expected, true, std::memory_order_acq_rel,
              std::memory_order_relaxed)) {
        next_bag_.insert(v);
      }
    }
  }

  void add_to_frontier(NodeId v) { add_to_frontier(v, load_dist(v)); }

  void visit_neighbors_parallel(NodeId u) {
    parallel_for(
        graph_.offsets[u], graph_.offsets[u + 1],
        [&](size_t i) {
          NodeId v = graph_.edges[i].v;
          NodeId dist_u = load_dist(u);
          NodeId new_dist = dist_u + 1;
          if (write_min_std(&dist_[v], new_dist)) {
            add_to_frontier(v, new_dist);
          }
        },
        BETA);
  }

  void visit_neighbors_sequential(NodeId u, NodeId *local_queue, size_t &rear) {
    for (EdgeId i = graph_.offsets[u]; i < graph_.offsets[u + 1]; i++) {
      NodeId v = graph_.edges[i].v;
      NodeId dist_u = load_dist(u);
      NodeId new_dist = dist_u + 1;
      if (write_min_std(&dist_[v], new_dist)) {
        if (rear < MAX_QUEUE_SIZE && new_dist < threshold_) {
          local_queue[rear++] = v;
        } else {
          add_to_frontier(v, new_dist);
        }
      }
    }
  }

  void sparse_relax(size_t frontier_size) {
    [[maybe_unused]] static const int num_threads = parlay::num_workers();
    const size_t queue_size = MAX_QUEUE_SIZE;
    parallel_for(0, frontier_size, [&](size_t i) {
      NodeId f = frontier_[i];
      assert(load_dist(f) < threshold_);
      if constexpr (use_local_queue) {
        NodeId local_queue[MAX_QUEUE_SIZE];
        size_t front = 0, rear = 0;
        size_t vertices_visited = 0;
        size_t edges_processed = 0;
        const size_t max_edges = queue_size * EDGE_PER_CACHELINE;
        local_queue[rear++] = f;
        while (front < rear && vertices_visited < MAX_QUEUE_SIZE &&
               edges_processed < max_edges) {
          NodeId u = local_queue[front++];
          assert(load_dist(u) < threshold_);
          vertices_visited++;
          size_t deg = graph_.offsets[u + 1] - graph_.offsets[u];
          edges_processed += deg;
          if (deg < BETA) {
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
      in_curr_frontier_[i].store(false, std::memory_order_relaxed);
      in_next_frontier_[i].store(false, std::memory_order_relaxed);
      dist_[i] = DIST_MAX;
    });
    dist_[s] = 0;
    in_curr_frontier_[s].store(true, std::memory_order_relaxed);
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
        parallel_for(0, frontier_size, [&](NodeId i) {
          in_curr_frontier_[frontier_[i]].store(false, std::memory_order_release);
        });
        sparse_relax(frontier_size);
        frontier_size = curr_bag_.pack_into(make_slice(frontier_));
      }
      round++;
      std::swap(curr_bag_, next_bag_);
      std::swap(in_curr_frontier_, in_next_frontier_);
      threshold_ += delta_;
    }
    return dist_;
  }
};
