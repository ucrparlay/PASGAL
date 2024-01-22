#pragma once
#include "parlay/monoid.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

template <class Seq, class Monoid>
class sparse_table {
  static constexpr int LOG2_BLOCK = 6;
  static constexpr size_t ST_BLOCK_SIZE = size_t{1} << LOG2_BLOCK;
  static constexpr size_t BLOCK_MASK = ST_BLOCK_SIZE - 1;

  const Seq &seq;
  Monoid m;
  using Type = decltype(m.identity);
  parlay::sequence<parlay::sequence<Type>> table;

 public:
  sparse_table(Seq &_seq, Monoid _m) : seq(_seq), m(_m) {
    size_t n = std::max(size_t{1}, seq.size() / ST_BLOCK_SIZE);
    size_t k = std::max(size_t{1}, (size_t)ceil(log2(n)));
    table =
        parlay::sequence<parlay::sequence<Type>>(k, parlay::sequence<Type>(n));
    parlay::parallel_for(0, n, [&](size_t i) {
      Type v = m.identity;
      for (size_t offset = 0;
           offset < ST_BLOCK_SIZE && ((i << LOG2_BLOCK) | offset) < seq.size();
           offset++) {
        v = m.f(v, seq[(i << LOG2_BLOCK) | offset]);
      }
      table[0][i] = v;
    });
    for (size_t i = 1; i < k; i++) {
      parlay::parallel_for(0, n, [&](size_t j) {
        if (j + (1 << i) <= n) {
          table[i][j] = m.f(table[i - 1][j], table[i - 1][j + (1 << (i - 1))]);
        }
      });
    }
  }
  Type query(size_t l, size_t r) {
    size_t block_l = (l >> LOG2_BLOCK) + 1, block_r = r >> LOG2_BLOCK;
    Type v = m.identity;
    if (block_l < block_r) {
      size_t s = 63 - __builtin_clzll(block_r - block_l);
      v = m.f(v, m.f(table[s][block_l], table[s][block_r - (1ull << s)]));
      for (size_t i = l; i < std::min(r, block_l << LOG2_BLOCK); i++) {
        v = m.f(v, seq[i]);
      }
      for (size_t i = std::max(l, block_r << LOG2_BLOCK); i < r; i++) {
        v = m.f(v, seq[i]);
      }
    } else {
      for (size_t i = l; i < r; i++) {
        v = m.f(v, seq[i]);
      }
    }
    return v;
  }
};