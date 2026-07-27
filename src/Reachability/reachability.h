#pragma once
#include <bit>
#include <climits>
#include <cstdint>

#include "chunkbag.h"
#include "graph.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "utils.h"

using namespace std;
using namespace parlay;

// Single-source reachability.  A vertex is either reachable or not, with no
// order to respect, so the local-queue walk costs no speculation and the pull
// pass stops at its first reachable in-neighbor.  `beta` bounds the walk's
// edges, trading rounds against redundant work.
template <class Graph>
class Reachability {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr size_t BLOCK_SIZE = 1024;
  static constexpr size_t SPARSE_TH = 20;
  static constexpr size_t MAX_QUEUE = 4096;
  static constexpr size_t BITS_PER_WORD = 64;

  const Graph &G;
  chunkbag<NodeId> bag;
  sequence<NodeId> frontier;
  sequence<uint8_t> visited;
  sequence<uint8_t> in_frontier;
  sequence<uint64_t> curr_bits;
  sequence<uint64_t> next_bits;
  size_t round_;
  size_t dense_rounds_;

  // A vertex is inserted only by the worker whose CAS marks it visited, so
  // it enters the bag at most once for the whole traversal.
  static size_t bag_capacity(size_t n) {
    return n + parlay::num_workers() * chunkbag<NodeId>::CHUNK_SIZE;
  }

  bool in_curr(NodeId v) const {
    return (curr_bits[v / BITS_PER_WORD] >> (v % BITS_PER_WORD)) & 1;
  }

 public:
  size_t beta = 2048;
  bool direction_optimizing = true;

  Reachability() = delete;
  Reachability(const Graph &_G) : G(_G), bag(bag_capacity(_G.n)) {
    frontier = sequence<NodeId>::uninitialized(G.n);
    visited = sequence<uint8_t>::uninitialized(G.n);
    in_frontier = sequence<uint8_t>::uninitialized(G.n);
    const size_t words = (G.n + BITS_PER_WORD - 1) / BITS_PER_WORD;
    curr_bits = sequence<uint64_t>(words, 0);
    next_bits = sequence<uint64_t>(words, 0);
  }

  size_t rounds() const { return round_; }
  size_t dense_rounds() const { return dense_rounds_; }

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
        if (rear < MAX_QUEUE) {
          local_queue[rear++] = v;
        } else {
          add_to_frontier(v);
        }
      }
    }
  }

  // Push: walk out-edges, extending locally until the edge budget is spent.
  void sparse_relax(size_t frontier_size) {
    parallel_for(
        0, frontier_size,
        [&](size_t i) {
          NodeId f = frontier[i];
          in_frontier[f] = false;
          NodeId local_queue[MAX_QUEUE];
          size_t front = 0, rear = 0, edges = 0;
          local_queue[rear++] = f;
          while (front < rear && edges < beta) {
            NodeId u = local_queue[front++];
            size_t deg = G.offsets[u + 1] - G.offsets[u];
            if (deg < BLOCK_SIZE) {
              visit_neighbors_sequential(u, local_queue, rear);
              edges += deg;
            } else {
              visit_neighbors_parallel(u);
            }
          }
          while (front < rear) {
            add_to_frontier(local_queue[front++]);
          }
        },
        1);
  }

  // Pull: every unvisited vertex looks for a frontier in-neighbor.  With no
  // distances to minimize, the first one found settles it, so the scan
  // always exits early.
  size_t dense_relax() {
    const size_t words = curr_bits.size();
    parallel_for(0, words, [&](size_t w) {
      uint64_t next = 0;
      const size_t base = w * BITS_PER_WORD;
      const size_t lim = std::min<size_t>(G.n - base, BITS_PER_WORD);
      for (size_t b = 0; b < lim; b++) {
        const NodeId u = (NodeId)(base + b);
        if (visited[u]) {
          continue;
        }
        const auto neighbors = G.in_neighors(u);
        for (size_t j = 0; j < neighbors.size(); j++) {
          if (in_curr(neighbors[j].v)) {
            visited[u] = true;
            next |= uint64_t(1) << b;
            break;
          }
        }
      }
      next_bits[w] = next;
    });
    std::swap(curr_bits, next_bits);
    return parlay::reduce(parlay::delayed_seq<size_t>(
        words, [&](size_t w) { return (size_t)std::popcount(curr_bits[w]); }));
  }

  // Move the packed frontier into the bitmap the pull pass reads.
  void sparse_to_dense(size_t frontier_size) {
    parallel_for(0, curr_bits.size(), [&](size_t w) { curr_bits[w] = 0; });
    parallel_for(0, frontier_size, [&](size_t i) {
      const NodeId v = frontier[i];
      in_frontier[v] = false;
      std::atomic_ref<uint64_t>(curr_bits[v / BITS_PER_WORD])
          .fetch_or(uint64_t(1) << (v % BITS_PER_WORD),
                    std::memory_order_relaxed);
    });
  }

  // Bitmap back to the bag, so the push direction can take over again.
  size_t dense_to_sparse() {
    parallel_for(0, curr_bits.size(), [&](size_t w) {
      uint64_t word = curr_bits[w];
      while (word) {
        const NodeId v = (NodeId)(w * BITS_PER_WORD + std::countr_zero(word));
        word &= word - 1;
        add_to_frontier(v);
      }
    });
    return bag.size();
  }

  sequence<uint8_t> reachability(NodeId s) {
    parallel_for(0, G.n, [&](size_t i) {
      in_frontier[i] = false;
      visited[i] = false;
    });
    parallel_for(0, curr_bits.size(), [&](size_t w) {
      curr_bits[w] = 0;
      next_bits[w] = 0;
    });
    bag.clear();

    round_ = 0;
    dense_rounds_ = 0;
    visited[s] = true;
    add_to_frontier(s);

    bool dense = false;
    size_t frontier_size = 1;
    while (frontier_size) {
      const bool want_dense =
          direction_optimizing && frontier_size * SPARSE_TH >= G.n;
      if (want_dense) {
        if (!dense) {
          frontier_size = bag.pack_into(make_slice(frontier));
          if (!frontier_size) {
            break;
          }
          sparse_to_dense(frontier_size);
          dense = true;
        }
        frontier_size = dense_relax();
        dense_rounds_++;
      } else {
        if (dense) {
          frontier_size = dense_to_sparse();
          dense = false;
          if (!frontier_size) {
            break;
          }
        }
        frontier_size = bag.pack_into(make_slice(frontier));
        if (!frontier_size) {
          break;
        }
        sparse_relax(frontier_size);
        frontier_size = bag.size();
      }
      round_++;
    }
    return visited;
  }
};
