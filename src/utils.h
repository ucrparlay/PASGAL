
#ifndef UTILS_H
#define UTILS_H

#include <atomic>
#include <functional>
#include <type_traits>

// Atomic primitives on plain objects, implemented with std::atomic_ref
// (C++20) instead of the previous __sync builtins.  Referenced objects
// must satisfy std::atomic_ref<ET>::required_alignment, which holds for
// the NodeId/EdgeId/bool sequences used across PASGAL.

template <typename ET>
inline bool atomic_compare_and_swap(ET *a, ET oldval, ET newval) {
  static_assert(sizeof(ET) <= 8, "Bad CAS length");
  return std::atomic_ref<ET>(*a).compare_exchange_strong(oldval, newval);
}

template <class ET>
inline bool compare_and_swap(std::atomic<ET> *a, ET oldval, ET newval) {
  return a->load(std::memory_order_relaxed) == oldval &&
         a->compare_exchange_strong(oldval, newval);
}

template <class ET>
inline bool compare_and_swap(ET *a, ET oldval, ET newval) {
  std::atomic_ref<ET> ref(*a);
  return ref.load(std::memory_order_relaxed) == oldval &&
         ref.compare_exchange_strong(oldval, newval);
}

template <typename E, typename EV>
inline E fetch_and_add(E *a, EV b) {
  if constexpr (std::is_integral_v<E>) {
    return std::atomic_ref<E>(*a).fetch_add(static_cast<E>(b));
  } else {
    std::atomic_ref<E> ref(*a);
    E oldV = ref.load(std::memory_order_relaxed);
    while (!ref.compare_exchange_weak(oldV, oldV + b));
    return oldV;
  }
}

template <typename E, typename EV>
inline void write_add(E *a, EV b) {
  (void)fetch_and_add(a, b);
}

template <typename ET, typename F = std::less<ET>>
inline bool write_min(ET *a, ET b, F less = {}) {
  std::atomic_ref<ET> ref(*a);
  ET c = ref.load(std::memory_order_relaxed);
  while (less(b, c)) {
    if (ref.compare_exchange_weak(c, b)) return true;
    // c holds the newly observed value after a failed CAS
  }
  return false;
}

// Like write_min, but on success reports the value it replaced in
// `old` (e.g. old == "infinity" identifies a first-time relaxation).
// `old` is unspecified when the function returns false.
template <typename ET, typename F = std::less<ET>>
inline bool write_min_old(ET *a, ET b, ET &old, F less = {}) {
  std::atomic_ref<ET> ref(*a);
  ET c = ref.load(std::memory_order_relaxed);
  while (less(b, c)) {
    if (ref.compare_exchange_weak(c, b)) {
      old = c;
      return true;
    }
  }
  return false;
}

template <typename ET, typename F = std::less<ET>>
inline bool write_max(ET *a, ET b, F less = {}) {
  std::atomic_ref<ET> ref(*a);
  ET c = ref.load(std::memory_order_relaxed);
  while (less(c, b)) {
    if (ref.compare_exchange_weak(c, b)) return true;
  }
  return false;
}

#endif  // UTILS_H
