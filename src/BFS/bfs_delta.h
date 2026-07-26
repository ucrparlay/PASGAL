#pragma once
#include <atomic>
#include <bit>
#include <climits>
#include <cstdint>
#include <functional>

#include "chunkbag.h"
#include "graph.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/utilities.h"

using namespace std;
using namespace parlay;

// Delta-stepping BFS with direction optimization.  Each band of `delta`
// levels repeats sub-iterations until empty, relaxing sparse (push) or
// dense (pull) by frontier density.  Sparse passes collect discoveries in
// a chunk bag; dense passes write a bitmap instead, which the pass can
// fill without atomics because it partitions the words.  Both carriers
// swap at the band boundary, so a dense band hands off to the next one
// without materializing a bag.
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
  // Direction switch: relax dense (pull) when the frontier holds at least
  // 1/SPARSE_THRESHOLD of all vertices, sparse (push) otherwise.  Same
  // single knob and value as bfs.h.
  static constexpr size_t SPARSE_THRESHOLD = 20;
  static constexpr size_t BITS_PER_WORD = 64;
  static constexpr size_t DENSE_BLOCK_WORDS = 8;  // pack blocking factor

  const Graph &graph_;
  NodeId threshold_;
  NodeId delta_;
  chunkbag<NodeId> curr_bag_;
  chunkbag<NodeId> next_bag_;
  sequence<NodeId> frontier_;
  sequence<NodeId> dist_;
  sequence<std::atomic<bool>> in_curr_frontier_;
  sequence<std::atomic<bool>> in_next_frontier_;
  // A frontier may live in the bag, in the bitmap, or -- in a band that
  // switches direction -- partly in each, so its size is the sum.
  sequence<uint64_t> dense_frontier_;       // this band
  sequence<uint64_t> dense_next_frontier_;  // next band
  sequence<size_t> dense_block_counts_;
  size_t dense_size_ = 0;       // popcount(dense_frontier_)
  size_t dense_next_size_ = 0;  // popcount(dense_next_frontier_)
  NodeId dense_min_ = 0;        // lower bound on the frontier's min distance
  constexpr static bool use_local_queue = true;

  NodeId load_dist(NodeId v,
                   std::memory_order order = std::memory_order_acquire) {
    return std::atomic_ref<NodeId>(dist_[v]).load(order);
  }

  void store_dist(NodeId v, NodeId value) {
    std::atomic_ref<NodeId>(dist_[v]).store(value, std::memory_order_release);
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

  // A vertex enters a bag only when its flag flips false->true, and flags
  // are cleared only by the retire pass that follows a pack (which empties
  // the bag), so a bag holds at most n -- plus one partial chunk each.
  static size_t bag_capacity(size_t n) {
    return n + parlay::num_workers() * chunkbag<NodeId>::CHUNK_SIZE;
  }

 public:
  BFS() = delete;
  BFS(const Graph &graph, NodeId delta)
      : graph_(graph),
        delta_(delta),
        curr_bag_(bag_capacity(graph.n)),
        next_bag_(bag_capacity(graph.n)) {
    frontier_ = sequence<NodeId>::uninitialized(graph.n);
    dist_ = sequence<NodeId>::uninitialized(graph.n);
    in_curr_frontier_ = sequence<std::atomic<bool>>::uninitialized(graph.n);
    in_next_frontier_ = sequence<std::atomic<bool>>::uninitialized(graph.n);
    const size_t num_words = (graph.n + BITS_PER_WORD - 1) / BITS_PER_WORD;
    dense_frontier_ = sequence<uint64_t>(num_words, 0);
    dense_next_frontier_ = sequence<uint64_t>(num_words, 0);
    dense_block_counts_ = sequence<size_t>::uninitialized(
        (num_words + DENSE_BLOCK_WORDS - 1) / DENSE_BLOCK_WORDS);
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

  // One sub-iteration.  Sizing needs no packing: the bag counts in
  // O(#workers) and the bitmap keeps a running popcount.
  bool relax() {
    size_t bag_size = curr_bag_.size();
    const size_t frontier_size = bag_size + dense_size_;
    if (frontier_size == 0) {
      return false;
    }

    if (frontier_size * SPARSE_THRESHOLD >= graph_.n) {
      // The pull pass reads only dist_ and rewrites every bitmap word, so
      // the incoming frontier is simply dropped: its members are expanded
      // implicitly when their neighbours pull from them.
      curr_bag_.clear();
      dense_relax();
    } else {
      size_t sparse_size;
      if (dense_size_ == 0) {
        sparse_size = curr_bag_.pack_into(make_slice(frontier_));
        retire_flags(sparse_size);
      } else {
        if (bag_size > 0) {
          // Fold rather than concatenate: bits dedup, so a vertex found
          // by both a sparse and a pull pass is relaxed once.
          fold_into_bitmap(curr_bag_.pack_into(make_slice(frontier_)));
        }
        // Packs into frontier_ directly; no flags to set or retire.
        sparse_size = dense_to_sparse();
        dense_size_ = 0;
      }
      sparse_relax(sparse_size);
    }
    return true;
  }

  // Must run before any relaxation can re-enqueue: clearing inside the
  // relaxation loop loses updates that land on a not-yet-cleared flag.
  void retire_flags(size_t frontier_size) {
    parallel_for(0, frontier_size, [&](size_t i) {
      in_curr_frontier_[frontier_[i]].store(false, std::memory_order_release);
    });
  }

  // Fold the already-packed bag into this band's bitmap.  Concurrent
  // writers can share a word, so the bit set is atomic.
  void fold_into_bitmap(size_t k) {
    parallel_for(0, k, [&](size_t i) {
      const NodeId v = frontier_[i];
      in_curr_frontier_[v].store(false, std::memory_order_relaxed);
      std::atomic_ref<uint64_t>(dense_frontier_[v / BITS_PER_WORD])
          .fetch_or(uint64_t(1) << (v % BITS_PER_WORD),
                    std::memory_order_relaxed);
    });
  }

  // Bottom-up sub-iteration: every vertex that can still improve scans
  // its in-neighbors for a parent below the threshold.  Words are owned
  // by one worker, so dist_ and both bitmaps are written without atomics.
  // Nothing can drop below dense_min_ + 1, so the scan stops there.
  void dense_relax() {
    const NodeId floor_dist = threshold_ - delta_;
    const NodeId cutoff = dense_min_ + 1;
    const size_t num_words = dense_frontier_.size();
    std::atomic<NodeId> frontier_min(DIST_MAX);
    parallel_for(0, num_words, [&](size_t w) {
      uint64_t next = 0;     // this band's frontier
      uint64_t promote = 0;  // vertices crossing into the next band
      NodeId local_min = DIST_MAX;
      const size_t base = w * BITS_PER_WORD;
      const size_t lim = std::min<size_t>(graph_.n - base, BITS_PER_WORD);
      for (size_t b = 0; b < lim; b++) {
        const NodeId u = (NodeId)(base + b);
        const NodeId du = load_dist(u);
        if (du < floor_dist) {  // settled in an earlier band
          continue;
        }
        if (du < threshold_ &&
            in_curr_frontier_[u].load(std::memory_order_relaxed)) {
          // Sparse-era flag surviving the direction switch: retire it.
          in_curr_frontier_[u].store(false, std::memory_order_relaxed);
        }
        if (du == floor_dist) {  // already minimal for this band
          continue;
        }
        NodeId new_du = du;
        const auto neighbors = graph_.in_neighors(u);
        for (size_t j = 0; j < neighbors.size(); j++) {
          NodeId dv = load_dist(neighbors[j].v);
          if (dv < threshold_ && dv + 1 < new_du) {
            new_du = dv + 1;
            if (new_du == cutoff) {
              break;
            }
          }
        }
        if (new_du != du) {
          store_dist(u, new_du);
          if (new_du < threshold_) {
            next |= uint64_t(1) << b;
            local_min = std::min(local_min, new_du);
          } else {
            promote |= uint64_t(1) << b;
          }
        }
      }
      dense_frontier_[w] = next;
      if (promote != 0) {  // accumulates across this band's passes
        dense_next_frontier_[w] |= promote;
      }
      // The minimum only increases and never falls below cutoff, so once
      // a word hits that floor the rest skip the CAS.
      if (local_min != DIST_MAX) {
        NodeId cur = frontier_min.load(std::memory_order_relaxed);
        while (local_min < cur &&
               !frontier_min.compare_exchange_weak(
                   cur, local_min, std::memory_order_relaxed)) {
        }
      }
    });
    dense_size_ = parlay::reduce(parlay::delayed_seq<size_t>(
        num_words,
        [&](size_t w) { return (size_t)std::popcount(dense_frontier_[w]); }));
    dense_next_size_ = parlay::reduce(parlay::delayed_seq<size_t>(
        num_words, [&](size_t w) {
          return (size_t)std::popcount(dense_next_frontier_[w]);
        }));
    const NodeId m = frontier_min.load(std::memory_order_relaxed);
    dense_min_ = (m == DIST_MAX) ? cutoff : m;
  }

  // Pack the bitmap straight into frontier_, clearing it on the way:
  // block popcounts scan into the output offsets.  No flag writes needed,
  // the pull passes already retired them.
  size_t dense_to_sparse() {
    const size_t num_words = dense_frontier_.size();
    const size_t num_blocks =
        (num_words + DENSE_BLOCK_WORDS - 1) / DENSE_BLOCK_WORDS;
    auto block_words = [&](size_t blk) {
      return std::make_pair(blk * DENSE_BLOCK_WORDS,
                            std::min(num_words, (blk + 1) * DENSE_BLOCK_WORDS));
    };
    parallel_for(0, num_blocks, [&](size_t blk) {
      auto [w_lo, w_hi] = block_words(blk);
      size_t count = 0;
      for (size_t w = w_lo; w < w_hi; w++) {
        count += (size_t)std::popcount(dense_frontier_[w]);
      }
      dense_block_counts_[blk] = count;
    });
    auto offsets = make_slice(dense_block_counts_).cut(0, num_blocks);
    const size_t total = parlay::scan_inplace(offsets);
    parallel_for(0, num_blocks, [&](size_t blk) {
      auto [w_lo, w_hi] = block_words(blk);
      size_t out = offsets[blk];
      for (size_t w = w_lo; w < w_hi; w++) {
        uint64_t word = dense_frontier_[w];
        if (word == 0) {
          continue;
        }
        dense_frontier_[w] = 0;
        do {
          frontier_[out++] =
              (NodeId)(w * BITS_PER_WORD + std::countr_zero(word));
          word &= word - 1;
        } while (word);
      }
    });
    return total;
  }

  void sparse_relax(size_t frontier_size) {
    parallel_for(0, frontier_size, [&](size_t i) {
      NodeId f = frontier_[i];
      assert(load_dist(f) < threshold_);
      if constexpr (use_local_queue) {
        // Walk several levels locally before the next global
        // sub-iteration -- what makes high-diameter graphs fast; shrinking
        // BETA to 128 edges costs a road graph 15%.  `rear` never exceeds
        // MAX_QUEUE_SIZE, so front < rear already caps the vertex count.
        NodeId local_queue[MAX_QUEUE_SIZE];
        size_t front = 0, rear = 0;
        size_t edges_processed = 0;
        local_queue[rear++] = f;
        while (front < rear && edges_processed < BETA) {
          NodeId u = local_queue[front++];
          assert(load_dist(u) < threshold_);
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

  // Split so a benchmark can time the solve without the reset.
  void prepare(NodeId s) {
    curr_bag_.clear();
    next_bag_.clear();
    parallel_for(0, graph_.n, [&](size_t i) {
      in_curr_frontier_[i].store(false, std::memory_order_relaxed);
      in_next_frontier_[i].store(false, std::memory_order_relaxed);
      dist_[i] = DIST_MAX;
    });
    parallel_for(0, dense_frontier_.size(), [&](size_t w) {
      dense_frontier_[w] = 0;
      dense_next_frontier_[w] = 0;
    });
    dense_size_ = dense_next_size_ = 0;
    dist_[s] = 0;
    in_curr_frontier_[s].store(true, std::memory_order_relaxed);
    threshold_ = delta_;
    dense_min_ = 0;
    curr_bag_.insert(s);
  }

  void run_prepared() {
    while (true) {
      bool band_had_work = false;
      while (relax()) {
        band_had_work = true;
      }
      // A discovery has dist = parent + 1 <= threshold_, so no vertex can
      // skip a band: an empty band means the traversal is done.
      if (!band_had_work) {
        break;
      }
      std::swap(curr_bag_, next_bag_);
      std::swap(in_curr_frontier_, in_next_frontier_);
      // Both carriers swap; the outgoing bitmap is all zeros because the
      // band ended on bag_size + dense_size_ == 0.
      std::swap(dense_frontier_, dense_next_frontier_);
      std::swap(dense_size_, dense_next_size_);
      threshold_ += delta_;
      dense_min_ = threshold_ - delta_;  // band floor: a valid lower bound
    }
  }

  const sequence<NodeId> &distances() const { return dist_; }

  sequence<NodeId> bfs(NodeId s) {
    prepare(s);
    run_prepared();
    return dist_;
  }
};
