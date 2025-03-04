
#ifndef SAMPLER_H
#define SAMPLER_H

#include <atomic>
#include <stdexcept>

using hash_t = uint32_t;

class Sampler {
  std::atomic<uint32_t> num_hits;
  uint32_t exp_hits;
  hash_t threshold;

 public:
  Sampler() {}

  Sampler(const size_t _exp_hits, double _sample_rate)
      : exp_hits(_exp_hits),
        threshold(_sample_rate * std::numeric_limits<hash_t>::max()) {
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

  void reset() { num_hits = 0; }

  uint32_t get_num_hits() const { return num_hits; }

  uint32_t get_exp_hits() const { return exp_hits; }

  void set_sample_rate(double _sample_rate) {
    threshold = _sample_rate * std::numeric_limits<hash_t>::max();
  }

  void reset(const size_t _exp_hits, double _sample_rate) {
    exp_hits = _exp_hits;
    threshold = _sample_rate * std::numeric_limits<hash_t>::max();
    num_hits = 0;
  }
};

#endif  // SAMPLER_H
