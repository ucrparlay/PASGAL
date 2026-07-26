
#ifndef CHUNKBAG_H
#define CHUNKBAG_H

#include <atomic>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

// Per-worker chunked bag: an unordered concurrent collection with the
// same interface as hashbag.  Each worker appends into a private chunk
// claimed from a shared pool (one fetch_add per CHUNK_SIZE insertions),
// so inserts are sequential cache-line writes with no hashing, no CAS
// retries, and no shared-counter traffic.  Every allocated chunk except
// each worker's current one is completely full (a worker claims a new
// chunk only after filling its old one), so pack copies dense chunks
// and consults only the <= num_workers() partial tails -- no empty-slot
// overscan and no clear pass.
template <class ET>
class chunkbag {
  static_assert(std::is_trivially_copyable_v<ET>,
                "chunkbag requires trivially copyable elements");

 public:
  static constexpr size_t CHUNK_SIZE = 1 << 11;  // elements per chunk

 private:
  struct alignas(64) WorkerState {
    size_t cur = 0, end = 0;  // unwritten range [cur, end); end==0: none
  };

  // The allocation cursor gets its own cache line: sharing one with
  // pool/states makes every chunk allocation invalidate the pointers that
  // every insert reads (profiled at ~12% of a low-diameter BFS).
  size_t capacity = 0;
  parlay::sequence<ET> pool;
  parlay::sequence<WorkerState> states;
  alignas(64) std::atomic<size_t> head{0};

 public:
  chunkbag() = default;

  // An exact slot count, not a hint: insert() aborts on overflow.  Only
  // the caller knows its insertion policy, so it owns the arithmetic --
  // usually (max live elements) + num_workers() * CHUNK_SIZE.
  explicit chunkbag(size_t _capacity)
      : capacity(_capacity),
        pool(parlay::sequence<ET>::uninitialized(capacity)),
        states(parlay::num_workers()),
        head(0) {}

  chunkbag(const chunkbag &other)
      : capacity(other.capacity),
        pool(other.pool),
        states(other.states),
        head(other.head.load()) {}

  chunkbag(chunkbag &&other) noexcept
      : capacity(other.capacity),
        pool(std::move(other.pool)),
        states(std::move(other.states)),
        head(other.head.load()) {}

  chunkbag &operator=(const chunkbag &other) {
    if (this != &other) {
      capacity = other.capacity;
      head.store(other.head.load());
      pool = other.pool;
      states = other.states;
    }
    return *this;
  }

  chunkbag &operator=(chunkbag &&other) noexcept {
    if (this != &other) {
      capacity = other.capacity;
      head.store(other.head.load());
      pool = std::move(other.pool);
      states = std::move(other.states);
    }
    return *this;
  }

  void clear() {
    head.store(0, std::memory_order_relaxed);
    for (auto &st : states) {
      st = WorkerState();
    }
  }

  // Exact element count since the last clear()/pack (callers de-duplicate
  // before inserting): allocated chunks minus each worker's unwritten
  // tail.  O(#workers); call only between insert phases.
  size_t size() const {
    size_t s = head.load(std::memory_order_relaxed);
    for (const WorkerState &st : states) {
      s -= st.end - st.cur;
    }
    return s;
  }

  void insert(ET u) {
    WorkerState &st = states[parlay::worker_id()];
    if (st.cur == st.end) {
      st.cur = head.fetch_add(CHUNK_SIZE, std::memory_order_relaxed);
      st.end = st.cur + CHUNK_SIZE;
      if (st.end > capacity) {  // asserts are off at -O3
        fprintf(stderr, "chunkbag: pool overflow (capacity %zu)\n", capacity);
        abort();
      }
    }
    pool[st.cur++] = u;
  }

  template <typename Seq>
  size_t pack_into(Seq &&out) {
    size_t len = head.load(std::memory_order_relaxed);
    size_t num_chunks = len / CHUNK_SIZE;
    if (num_chunks == 0) {
      clear();
      return 0;
    }
    // Chunk fill counts: CHUNK_SIZE unless it is a worker's current
    // chunk, in which case only [chunk_begin, cur) is valid.
    auto fill = parlay::sequence<size_t>(num_chunks, CHUNK_SIZE);
    for (const WorkerState &st : states) {
      if (st.end != 0) {
        size_t chunk_begin = st.end - CHUNK_SIZE;
        fill[chunk_begin / CHUNK_SIZE] = st.cur - chunk_begin;
      }
    }
    auto dest = fill;  // copy; scan gives per-chunk output offsets
    size_t total = parlay::scan_inplace(dest);
    parlay::parallel_for(
        0, num_chunks,
        [&](size_t c) {
          size_t src = c * CHUNK_SIZE;
          size_t dst = dest[c];
          size_t cnt = fill[c];
          for (size_t i = 0; i < cnt; i++) {
            out[dst + i] = pool[src + i];
          }
        },
        1);
    clear();
    return total;
  }

  template <typename Seq, typename UnaryPred>
  size_t pack_into_pred(Seq &&out, UnaryPred &&f) {
    size_t len = head.load(std::memory_order_relaxed);
    size_t num_chunks = len / CHUNK_SIZE;
    if (num_chunks == 0) {
      clear();
      return 0;
    }
    auto fill = parlay::sequence<size_t>(num_chunks, CHUNK_SIZE);
    for (const WorkerState &st : states) {
      if (st.end != 0) {
        size_t chunk_begin = st.end - CHUNK_SIZE;
        fill[chunk_begin / CHUNK_SIZE] = st.cur - chunk_begin;
      }
    }
    auto counts = parlay::sequence<size_t>::uninitialized(num_chunks);
    parlay::parallel_for(
        0, num_chunks,
        [&](size_t c) {
          size_t src = c * CHUNK_SIZE, cnt = 0;
          for (size_t i = 0; i < fill[c]; i++) {
            if (f(pool[src + i])) cnt++;
          }
          counts[c] = cnt;
        },
        1);
    size_t total = parlay::scan_inplace(counts);
    parlay::parallel_for(
        0, num_chunks,
        [&](size_t c) {
          size_t src = c * CHUNK_SIZE;
          size_t o = counts[c];
          for (size_t i = 0; i < fill[c]; i++) {
            if (f(pool[src + i])) out[o++] = pool[src + i];
          }
        },
        1);
    clear();
    return total;
  }
};

#endif  // CHUNKBAG_H
