// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010-2016 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include <tuple>
#include <unordered_set>
#include <vector>

#include "parlay/primitives.h"
#include "parlay/sequence.h"

using namespace parlay;

namespace gbbs {
inline size_t hashToRange(const size_t& h, const size_t& mask) {
  return h & mask;
}
inline size_t incrementIndex(const size_t& h, const size_t& mask) {
  return hashToRange(h + 1, mask);
}

template <class K, class KeyHash>
class resizable_table {
 public:
  size_t m;
  size_t mask;
  K empty;
  sequence<K> table;
  KeyHash key_hash;

  inline size_t firstIndex(K& k) { return hashToRange(key_hash(k), mask); }

  resizable_table() : m(0) { mask = 0; }

  resizable_table(size_t _m, K _empty, KeyHash _key_hash)
      : m((size_t)1 << parlay::log2_up((size_t)(_m * 1.2))),
        mask(m - 1),
        empty(_empty),
        key_hash(_key_hash) {
    table = sequence<K>(m, empty);
  }

  void insert(K k) {
    size_t h = firstIndex(k);
    while (true) {
      if (table[h] == empty && atomic_compare_and_swap(&table[h], empty, k)) {
        return;
      }
      h = incrementIndex(h, mask);
    }
  }

  sequence<K> entries() {
    auto pred = [&](K& t) { return t != empty; };
    return filter(make_slice(table), pred);
  }
};

template <class K, class KeyHash>
inline resizable_table<K, KeyHash> make_resizable_table(size_t m, K empty,
                                                        KeyHash key_hash) {
  return resizable_table<K, KeyHash>(m, empty, key_hash);
}

}  // namespace gbbs