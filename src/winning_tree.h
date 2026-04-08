#ifndef WINNING_TREE_H
#define WINNING_TREE_H

#include <array>
#include <bit>
#include <atomic>
#include <cassert>
#include <cstdint>

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

// A 64-ary hierarchical bitmap for efficient concurrent insertion and iteration
// of set bit indices. Each level uses atomic uint64_t words where each bit
// represents the bitwise OR of 64 bits from the level below.
//
// Performance optimizations:
// - Uses memory_order_relaxed for atomic operations (phase-concurrent)
// - Early termination when bit is already set
// - Cache-friendly iteration with tree pruning
// - C++20 standard bit operations
class WinningTree {
 private:
  static constexpr size_t BITS_PER_WORD = 64;

  size_t n_;                                      // Total number of bits
  size_t num_levels_;                             // Height of the tree
  parlay::sequence<std::atomic<uint64_t>> data_;  // Contiguous array for all levels
  parlay::sequence<size_t> level_offsets_;        // Starting offset of each level
  parlay::sequence<size_t> level_sizes_;          // Number of words per level

 public:
  // Constructs a winning tree that can represent n bits
  explicit WinningTree(size_t n) : n_(n) {
    assert(n > 0);

    // Calculate the size of each level bottom-up
    // Level 0 (bottom): ceil(n / 64) words
    // Level i+1: ceil(level[i].size / 64) words
    parlay::sequence<size_t> sizes;
    size_t curr = (n + BITS_PER_WORD - 1) / BITS_PER_WORD;
    sizes.push_back(curr);
    while (curr > 1) {
      curr = (curr + BITS_PER_WORD - 1) / BITS_PER_WORD;
      sizes.push_back(curr);
    }

    num_levels_ = sizes.size();
    level_sizes_ = sizes;

    // Calculate offsets and total size
    level_offsets_ = parlay::sequence<size_t>(num_levels_);
    size_t total_words = 0;
    for (size_t i = 0; i < num_levels_; i++) {
      level_offsets_[i] = total_words;
      total_words += level_sizes_[i];
    }

    // Allocate all levels contiguously for better cache locality
    data_ = parlay::sequence<std::atomic<uint64_t>>(total_words);
    for (size_t i = 0; i < total_words; i++) {
      data_[i].store(0, std::memory_order_relaxed);
    }
  }

  // Non-copyable for safety with atomic types
  WinningTree(const WinningTree&) = delete;
  WinningTree& operator=(const WinningTree&) = delete;
  WinningTree(WinningTree&&) noexcept = default;
  WinningTree& operator=(WinningTree&&) noexcept = default;

  // Sets the bit at the given index. Thread-safe for concurrent inserts.
  // Propagates the change up the tree using atomic OR operations.
  void insert(size_t index) {
    assert(index < n_);

    size_t pos = index;
    for (size_t level = 0; level < num_levels_; level++) {
      size_t word_idx = pos / BITS_PER_WORD;
      size_t bit_idx = pos % BITS_PER_WORD;
      uint64_t mask = 1ULL << bit_idx;

      // Atomically set the bit and get the old value
      uint64_t old_val = data_[level_offsets_[level] + word_idx].fetch_or(
          mask, std::memory_order_relaxed);

      // Early exit: if bit was already set, no need to propagate further
      if (old_val & mask) {
        return;
      }

      // Move to parent in next level
      pos = word_idx;
    }
  }

  bool empty() const {
    return data_[level_offsets_[num_levels_ - 1]].load(
               std::memory_order_relaxed) == 0;
  }

  void clear() {
    parlay::parallel_for(0, data_.size(), [&](size_t i) {
      data_[i].store(0, std::memory_order_relaxed);
    });
  }

  // Parallel top-down iteration over all set bit indices in [0, n).
  // Notes:
  // 1) func(index) must be thread-safe.
  // 2) Visit order is non-deterministic.
  // 3) Sparse workloads benefit from pruning empty subtrees before reaching leaves.
  template <typename Func>
  void iterate_all(Func&& func) {
    if (empty()) {
      return;
    }
    auto&& fn = func;
    iterate_level_parallel<false>(num_levels_ - 1, 0, fn);
  }

  // Iterates over a snapshot of the current tree and clears visited nodes as it goes.
  template <typename Func>
  void iterate_all_and_clear(Func&& func) {
    if (empty()) {
      return;
    }
    auto&& fn = func;
    iterate_level_parallel<true>(num_levels_ - 1, 0, fn);
  }

  // Returns the number of set bits.
  // Not thread-safe with concurrent inserts (provides approximate count).
  size_t size() const {
    return parlay::reduce(parlay::tabulate<size_t>(level_sizes_[0], [&](size_t i) {
      uint64_t word = data_[level_offsets_[0] + i].load(std::memory_order_relaxed);
      return static_cast<size_t>(std::popcount(word));
    }));
  }

 private:
  // Iterate over set bits in a word, calling action(base + bit_position)
  // for each set bit.
  template <typename Action>
  static void for_each_bit(uint64_t word, size_t base, Action&& action) {
    while (word) {
      size_t bit = std::countr_zero(word);
      action(base + bit);
      word &= word - 1;
    }
  }

  template <bool Clear, typename Func>
  void iterate_level_parallel(size_t level, size_t word_idx, Func& func) {
    uint64_t word;
    if constexpr (Clear) {
      word = data_[level_offsets_[level] + word_idx].exchange(
          0, std::memory_order_relaxed);
    } else {
      word = data_[level_offsets_[level] + word_idx].load(std::memory_order_relaxed);
    }

    if (word == 0) {
      return;
    }

    const size_t base = word_idx << 6;

    if (level == 0) {
      for_each_bit(word, base, func);
      return;
    }

    std::array<size_t, BITS_PER_WORD> child_indices;
    size_t child_count = 0;
    for_each_bit(word, base, [&](size_t idx) {
      child_indices[child_count++] = idx;
    });

    // At the topmost level, subtree sizes can be heavily skewed, so use
    // granularity=1 to fork every child and let work-stealing balance.
    // At lower levels, subtrees are bounded and more uniform, so let
    // parlay's granularity speculation (granularity=0) decide.
    parlay::parallel_for(
        0, child_count,
        [&](size_t i) {
          iterate_level_parallel<Clear>(level - 1, child_indices[i], func);
        },
        level == num_levels_ - 1 ? 1 : 0);
  }

};

#endif  // WINNING_TREE_H
