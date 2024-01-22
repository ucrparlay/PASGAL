
#ifndef SAMPLER_H
#define SAMPLER_H

#include <atomic>
using hash_t = uint64_t;

class Sampler {
  std::atomic<uint32_t> num_hits;
  const uint32_t exp_hits;
  hash_t threshold;

 public:
  Sampler(const size_t _exp_hits, hash_t _threshold)
      : exp_hits(_exp_hits), threshold(_threshold) {
    num_hits = 0;
  }

  Sampler(const Sampler &other)
      : exp_hits(other.exp_hits), threshold(other.threshold) {
    num_hits = other.num_hits.load();
  }

  Sampler(Sampler &&other)
      : exp_hits(other.exp_hits), threshold(other.threshold) {
    num_hits = other.num_hits.load();
  }

  bool sample(hash_t random_number, bool &callback) {
    callback = false;
    if (num_hits >= exp_hits) {
      return false;
    }
    if (random_number < threshold) {
      uint32_t ret = num_hits.fetch_add(1);
      if (ret >= exp_hits) {
        return false;
      } else if (ret + 1 == exp_hits) {
        callback = true;
      }
    }
    return true;
  }

  void set_thres(hash_t _threshold) { threshold = _threshold; }

  void reset() { num_hits = 0; }
};

#endif  // SAMPLER_H
