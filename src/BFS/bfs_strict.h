#pragma once

#include <atomic>
#include <bit>
#include <cstdint>
#include <limits>

#include "chunkbag.h"
#include "graph.h"
#include "hashbag.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"

using namespace std;
using namespace parlay;

// Strict level-synchronous BFS with a chunked sparse frontier and a bitmap
// pull frontier. The direction policy matches GBBS: use pull when
// |F| + edges(F) > m/20, and keep a dense frontier dense while |F| > n/10.
template <class Graph, template <class> class BagT = chunkbag>
class BFS {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr NodeId DIST_MAX = numeric_limits<NodeId>::max();
  static constexpr size_t EDGE_PARALLEL_THRESHOLD = 2048;
  static constexpr size_t EDGE_GRAIN = 1024;
  static constexpr size_t DENSE_ALPHA = 20;
  static constexpr size_t DENSE_STAY_BETA = 10;
  static constexpr size_t BLOCK_SIZE = 2048;
  static constexpr size_t BITS_PER_WORD = 64;
  static_assert(BLOCK_SIZE % BITS_PER_WORD == 0);

  struct alignas(64) WorkerEdgeCount {
    size_t value = 0;
  };

  const Graph &graph_;
  BagT<NodeId> next_bag_;
  sequence<NodeId> frontier_;
  sequence<NodeId> dist_;
  sequence<uint64_t> dense_frontier_;
  sequence<uint64_t> dense_next_frontier_;
  sequence<size_t> block_counts_;
  sequence<size_t> block_offsets_;
  sequence<WorkerEdgeCount> sparse_edge_counts_;
  size_t frontier_size_ = 0;
  size_t frontier_edges_ = 0;
  NodeId level_ = 0;
  bool dense_ = false;

  bool discover(NodeId v, NodeId distance) {
    atomic_ref<NodeId> ref(dist_[v]);
    if (ref.load(memory_order_relaxed) != DIST_MAX) {
      return false;
    }
    NodeId expected = DIST_MAX;
    return ref.compare_exchange_strong(expected, distance,
                                       memory_order_relaxed,
                                       memory_order_relaxed);
  }

  void visit_edge(EdgeId edge, NodeId next_distance) {
    const NodeId v = graph_.edges[edge].v;
    if (discover(v, next_distance)) {
      next_bag_.insert(v);
      sparse_edge_counts_[worker_id()].value +=
          static_cast<size_t>(graph_.offsets[v + 1] -
                              graph_.offsets[v]);
    }
  }

  void visit_neighbors(NodeId u, NodeId next_distance) {
    const EdgeId begin = graph_.offsets[u];
    const EdgeId end = graph_.offsets[u + 1];
    if (end - begin < EDGE_PARALLEL_THRESHOLD) {
      for (EdgeId edge = begin; edge < end; edge++) {
        visit_edge(edge, next_distance);
      }
      return;
    }
    parallel_for(
        begin, end,
        [&](EdgeId edge) { visit_edge(edge, next_distance); }, EDGE_GRAIN);
  }

  size_t sparse_relax(size_t frontier_size, NodeId next_distance) {
    for (auto &count : sparse_edge_counts_) {
      count.value = 0;
    }
    parallel_for(
        0, frontier_size,
        [&](size_t i) { visit_neighbors(frontier_[i], next_distance); }, 1);
    const size_t next_size = next_bag_.pack_into(make_slice(frontier_));
    size_t next_edges = 0;
    for (const auto &count : sparse_edge_counts_) {
      next_edges += count.value;
    }
    frontier_edges_ = next_edges;
    return next_size;
  }

  bool should_use_dense(size_t frontier_size) const {
    return frontier_size + frontier_edges_ > graph_.m / DENSE_ALPHA;
  }

  void build_dense_frontier(size_t frontier_size) {
    parallel_for(0, dense_frontier_.size(),
                 [&](size_t word) { dense_frontier_[word] = 0; }, 1);
    parallel_for(
        0, frontier_size,
        [&](size_t i) {
          const size_t v = frontier_[i];
          const size_t word = v / BITS_PER_WORD;
          const uint64_t bit = uint64_t{1} << (v % BITS_PER_WORD);
          atomic_ref<uint64_t>(dense_frontier_[word])
              .fetch_or(bit, memory_order_relaxed);
        },
        1);
  }

  bool in_dense_frontier(NodeId v) const {
    const size_t index = static_cast<size_t>(v);
    return (dense_frontier_[index / BITS_PER_WORD] >>
            (index % BITS_PER_WORD)) &
           uint64_t{1};
  }

  size_t dense_relax(NodeId next_distance) {
    const size_t num_blocks = (graph_.n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    parallel_for(
        0, num_blocks,
        [&](size_t block) {
          const size_t begin = block * BLOCK_SIZE;
          const size_t end = min(graph_.n, begin + BLOCK_SIZE);
          const size_t first_word = begin / BITS_PER_WORD;
          const size_t last_word =
              (end + BITS_PER_WORD - 1) / BITS_PER_WORD;
          size_t count = 0;
          for (size_t word_index = first_word; word_index < last_word;
               word_index++) {
            const size_t word_begin =
                max(begin, word_index * BITS_PER_WORD);
            const size_t word_end =
                min(end, (word_index + 1) * BITS_PER_WORD);
            uint64_t next_word = 0;
            for (size_t index = word_begin; index < word_end; index++) {
              if (dist_[index] != DIST_MAX) {
                continue;
              }
              const auto neighbors =
                  graph_.in_neighors(static_cast<NodeId>(index));
              for (size_t j = 0; j < neighbors.size(); j++) {
                if (in_dense_frontier(neighbors[j].v)) {
                  dist_[index] = next_distance;
                  next_word |= uint64_t{1} << (index % BITS_PER_WORD);
                  break;
                }
              }
            }
            dense_next_frontier_[word_index] = next_word;
            count += popcount(next_word);
          }
          block_counts_[block] = count;
        },
        1);
    size_t next_size = 0;
    for (size_t block = 0; block < num_blocks; block++) {
      next_size += block_counts_[block];
    }
    swap(dense_frontier_, dense_next_frontier_);
    return next_size;
  }

  size_t materialize_dense_frontier() {
    const size_t num_blocks = (graph_.n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    parallel_for(
        0, num_blocks,
        [&](size_t block) {
          const size_t begin = block * BLOCK_SIZE;
          const size_t end = min(graph_.n, begin + BLOCK_SIZE);
          size_t count = 0;
          for (size_t index = begin; index < end; index += BITS_PER_WORD) {
            count += popcount(dense_frontier_[index / BITS_PER_WORD]);
          }
          block_offsets_[block] = count;
        },
        1);
    auto offsets =
        make_slice(block_offsets_.begin(), block_offsets_.begin() + num_blocks);
    const size_t total = scan_inplace(offsets);
    parallel_for(
        0, num_blocks,
        [&](size_t block) {
          const size_t begin = block * BLOCK_SIZE;
          const size_t end = min(graph_.n, begin + BLOCK_SIZE);
          size_t output = block_offsets_[block];
          size_t edge_count = 0;
          for (size_t index = begin; index < end; index += BITS_PER_WORD) {
            uint64_t word = dense_frontier_[index / BITS_PER_WORD];
            while (word != 0) {
              const size_t bit = countr_zero(word);
              const size_t v = index + bit;
              if (v < graph_.n) {
                frontier_[output++] = static_cast<NodeId>(v);
                edge_count +=
                    static_cast<size_t>(graph_.offsets[v + 1] -
                                        graph_.offsets[v]);
              }
              word &= word - 1;
            }
          }
          block_counts_[block] = edge_count;
        },
        1);
    size_t edge_count = 0;
    for (size_t block = 0; block < num_blocks; block++) {
      edge_count += block_counts_[block];
    }
    frontier_edges_ = edge_count;
    return total;
  }

 public:
  BFS() = delete;

  // Bag sizing lives here, not in the bag: a vertex is inserted only by
  // the worker whose CAS moves its distance off DIST_MAX, which succeeds
  // exactly once per vertex for the whole traversal, so at most n
  // elements are ever live.  Chunk bags leave each worker one partial
  // chunk on top of that.
  static size_t bag_capacity(size_t n) {
    if constexpr (requires { BagT<NodeId>::CHUNK_SIZE; }) {
      return n + parlay::num_workers() * BagT<NodeId>::CHUNK_SIZE;
    } else {
      return n;  // hashbag sizes its own table from the element count
    }
  }

  explicit BFS(const Graph &graph)
      : graph_(graph),
        next_bag_(bag_capacity(graph.n)),
        frontier_(sequence<NodeId>::uninitialized(graph.n)),
        dist_(sequence<NodeId>::uninitialized(graph.n)),
        dense_frontier_((graph.n + BITS_PER_WORD - 1) / BITS_PER_WORD, 0),
        dense_next_frontier_((graph.n + BITS_PER_WORD - 1) / BITS_PER_WORD, 0),
        block_counts_((graph.n + BLOCK_SIZE - 1) / BLOCK_SIZE),
        block_offsets_((graph.n + BLOCK_SIZE - 1) / BLOCK_SIZE),
        sparse_edge_counts_(num_workers()) {}

  void prepare(NodeId source) {
    next_bag_.clear();
    parallel_for(0, graph_.n, [&](size_t i) { dist_[i] = DIST_MAX; });
    dist_[source] = 0;
    frontier_[0] = source;
    frontier_size_ = 1;
    frontier_edges_ =
        static_cast<size_t>(graph_.offsets[source + 1] -
                            graph_.offsets[source]);
    level_ = 0;
    dense_ = false;
  }

  void run_prepared() {
    while (frontier_size_ != 0) {
      const NodeId next_distance = level_ + 1;
      bool use_dense =
          dense_ && frontier_size_ > graph_.n / DENSE_STAY_BETA;
      const bool bitmap_valid = dense_;

      if (!use_dense) {
        if (dense_) {
          frontier_size_ = materialize_dense_frontier();
          dense_ = false;
        }
        use_dense = should_use_dense(frontier_size_);
      }

      if (use_dense) {
        if (!bitmap_valid) {
          build_dense_frontier(frontier_size_);
        }
        dense_ = true;
        frontier_size_ = dense_relax(next_distance);
      } else {
        frontier_size_ = sparse_relax(frontier_size_, next_distance);
      }
      level_ = next_distance;
    }
  }

  const sequence<NodeId> &distances() const { return dist_; }
};
