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
#include <cstddef>

#include <algorithm>
#include <atomic>
#include <optional>
#include <utility>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

template <typename K>
struct iter_k{
    using index = unsigned long;
    K k;
    index i;
    size_t num_prob;
};

template <typename K,typename V,typename Hash = parlay::hash<K>>
struct hash_table {
 private:
  static constexpr size_t kResizableTableCacheLineSz = 128;

 public:
  using KV = std::pair<K,V>;
  using index = unsigned long;
  size_t m;
  size_t mask;
  KV empty;
  index first_index(K k) { return hash(k) & mask;}
  index next_index(index h) { return (h + 1)& mask;}
  Hash hash;
  bool equal(KV a, KV b){return 
    parlay::internal::get_key(a)==parlay::internal::get_key(b) &&
    parlay::internal::get_val(a)==parlay::internal::get_val(b);}
  parlay::sequence<KV> H;
  parlay::sequence<size_t> cts;
  size_t ne;
  bool overfull;

  hash_table(long size,  KV _empty, Hash&& _hash = {}) 
    // : m((size_t)1 << parlay::log2_up((size_t)(2*size))),
    : m((size_t)1 << parlay::log2_up((size_t)(size))),
      mask(m-1),empty(_empty),hash(_hash),ne(0), overfull(false)
      {
      H = sequence<KV>::uninitialized(m);
      size_t workers = num_workers();
      cts = sequence<size_t>::uninitialized(kResizableTableCacheLineSz * workers);
      parallel_for(0, m, [&](size_t i){H[i]=empty;});
      for (size_t i = 0; i < workers; i++) {
        cts[i * kResizableTableCacheLineSz] = 0;
      } 
  }
  void double_size() {
    m = 4 * m;
    mask = m - 1;
    size();
    ne = 0;
    H = sequence<KV>::uninitialized(m);
    parallel_for(0, m, [&](size_t i){H[i]=empty;});
    overfull = false;
  }

  size_t size(){
    for (size_t i = 0;i<parlay::num_workers(); i++){
        ne+= cts[i*kResizableTableCacheLineSz];
        cts[i*kResizableTableCacheLineSz]=0;
    }
    if (ne >= m){overfull=true;}
    return ne;
  }

  bool insert(const K& k, const V& v) {
    KV kv = std::make_pair(k,v);
    index i = first_index(k);
    for (size_t count=0; count < 2000; count++) {
      if (equal(H[i], empty) && atomic_compare_and_swap(&H[i],empty,kv)){
        size_t wn = parlay::worker_id();
        cts[wn * kResizableTableCacheLineSz]++;
        return true;
      }
      if (equal(H[i], kv)) return false;
	    i = next_index(i);
	  }
    // std::cout << "Hash table overfull" << std::endl;
    overfull=true;
    return false;
  }

  bool contain(const KV kv) {
    K k = parlay::internal::get_key(kv);
    index i = first_index(k);
    while (true) {
      if (equal(H[i], empty)) return false;
      if (equal(H[i], kv)) return true;
      i = next_index(i);
    }
  }

  template <class F>
  void map(F& f) {
    parlay::parallel_for(0, m, [&](size_t i) {
      if (!equal(H[i], empty)) {
        f(H[i]);
      }
    });
  }

  bool init_iter(iter_k<K>& iter){
    index i = first_index(iter.k);
    while (true){
      if (equal(H[i],empty)){
        return false;
      }if (parlay::internal::get_key(H[i]) == iter.k){
        iter.num_prob=0;
        iter.i=i;
        return true;
      }
      i = next_index(i);
    }
  }
  bool has_next(iter_k<K>& iter){
    while (iter.num_prob<m){
      iter.i = next_index(iter.i);
      iter.num_prob++;
      if (equal(H[iter.i],empty)){
        return false;
      }if (parlay::internal::get_key(H[iter.i]) ==iter.k){
        return true;
      }
    }
    return false;
  }
  V get(iter_k<K>& iter){
    return parlay::internal::get_val(H[iter.i]);
  }
  parlay::sequence<KV> pack(){
    return parlay::filter(H, [&](auto p){return p != empty;});
  }
  // parlay::sequence<V> get_values(K key){
  //   iter_k<K> u_iter{key, 0,0};
  //   parlay::sequence<V> values;
  //   if (init_iter(u_iter)){
  //     while(true){
  //       V u_label=get(u_iter);
  //       values.push_back(u_label);
  //       if (!has_next(u_iter)){break;}
  //     }
  //   }
  //   return values;
  // }
};