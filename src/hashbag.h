
#ifndef HASHBAG_H
#define HASHBAG_H

#include <atomic>

#include "parlay/delayed_sequence.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"
#include "sampler.h"
#include "utils.h"

template <class ET>
class hashbag {
  static constexpr size_t BLOCK_SIZE = 1 << 10;
  static constexpr size_t MIN_BAG_SIZE = 1 << 6;
  static constexpr size_t MAX_PROBES = 1000;
  static constexpr size_t EXP_NUM_SAMPLES = 32;

  size_t n;
  const ET empty;
  std::atomic<uint32_t> bag_id;

  parlay::sequence<size_t> bag_sizes;
  parlay::sequence<size_t> offsets;
  parlay::sequence<Sampler> samplers;
  parlay::sequence<ET> pool;

 public:
  hashbag() = default;

  hashbag(size_t _n, double load_factor = 0.5,
          const ET _empty = std::numeric_limits<ET>::max())
      : n(_n), empty(_empty) {
    bag_id = 0;
    size_t cur_size = MIN_BAG_SIZE;
    size_t total_size = 0;
    // for (size_t i = 0; total_size * load_factor < n; i++) {
    while (total_size * load_factor < n) {
      double sample_rate = EXP_NUM_SAMPLES / (cur_size * load_factor);
      bag_sizes.push_back(cur_size);
      offsets.push_back(total_size);
      samplers.push_back(Sampler(EXP_NUM_SAMPLES, sample_rate));
      total_size += cur_size;
      cur_size *= 2;
    }
    pool = parlay::sequence<ET>(total_size, empty);
  }

  hashbag(const hashbag &other)
      : n(other.n),
        empty(other.empty),
        bag_id(other.bag_id.load()),
        bag_sizes(other.bag_sizes),
        offsets(other.offsets),
        samplers(other.samplers),
        pool(other.pool) {}

  hashbag(hashbag &&other)
      : n(other.n),
        empty(other.empty),
        bag_id(other.bag_id.load()),
        bag_sizes(other.bag_sizes),
        offsets(other.offsets),
        samplers(other.samplers),
        pool(other.pool) {}

  void clear() {
    for (size_t i = 0; i <= bag_id; i++) {
      samplers[i].reset();
    }
    parlay::parallel_for(
        0, offsets[bag_id] + bag_sizes[bag_id],
        [&](size_t i) { pool[i] = empty; }, BLOCK_SIZE);
    bag_id = 0;
  }

  void insert(ET u) {
    uint32_t local_id = bag_id;
    auto random_number = parlay::hash32(u);
    size_t idx = random_number & (bag_sizes[local_id] - 1);
    bool callback = false;
    while (local_id + 1 < bag_sizes.size() &&
           !samplers[local_id].sample(random_number, callback)) {
      compare_and_swap(&bag_id, local_id, local_id + 1);
      local_id = bag_id;
    }
    size_t num_probes = 0;
    idx = random_number & (bag_sizes[local_id] - 1);
    while (!compare_and_swap(&pool[offsets[local_id] + idx], empty, u)) {
      idx++;
      if (idx == bag_sizes[local_id]) {
        idx = 0;
      }
      num_probes++;
      if (num_probes == bag_sizes[local_id] || num_probes == MAX_PROBES) {
        num_probes = 0;
        compare_and_swap(&bag_id, local_id, local_id + 1);
        if (local_id != bag_id) {
          local_id = bag_id;
          assert(local_id < bag_sizes.size() && "hashbag is full");
          idx = random_number & (bag_sizes[local_id] - 1);
        }
      }
    }
  }

  parlay::sequence<ET> pack() {
    size_t len = offsets[bag_id] + bag_sizes[bag_id];
    auto pred = parlay::delayed_seq<bool>(
        len, [&](size_t i) { return pool[i] != empty; });
    parlay::sequence<ET> records = parlay::pack(pool.cut(0, len), pred);
    clear();
    return records;
  }

  template <typename Seq>
  size_t pack_into(Seq &&out) {
    size_t len = offsets[bag_id] + bag_sizes[bag_id];
    auto pred = parlay::delayed_seq<bool>(
        len, [&](size_t i) { return pool[i] != empty; });
    size_t num_records = parlay::pack_into_uninitialized(pool.cut(0, len), pred,
                                                         make_slice(out));
    clear();
    return num_records;
  }

  template <typename Seq, typename UnaryPred>
  size_t pack_into_pred(Seq &&out, UnaryPred &&f) {
    size_t len = offsets[bag_id] + bag_sizes[bag_id];
    auto pred = parlay::delayed_seq<bool>(
        len, [&](size_t i) { return pool[i] != empty && f(pool[i]); });
    size_t num_records = parlay::pack_into_uninitialized(pool.cut(0, len), pred,
                                                         make_slice(out));
    clear();
    return num_records;
  }

  void print() {
    for (size_t i = 0; i <= bag_id; i++) {
      auto seq = parlay::delayed_seq<uint32_t>(bag_sizes[i], [&](size_t j) {
        return pool[j + offsets[i]] != empty;
      });
      size_t ret = parlay::reduce(seq);
      printf("i=%lu: size=%zu, capacity=%zu, load=%f\n", i, ret, bag_sizes[i],
             1.0 * ret / bag_sizes[i]);
    }
  }
};

#endif  // HASHBAG_H
