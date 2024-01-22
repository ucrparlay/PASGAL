#pragma once

/* Union-Find options */
enum FindOption {
  find_compress,
  find_naive,
  find_split,
  find_halve,
  find_atomic_split,
  find_atomic_halve
};
enum UniteOption { unite, unite_early, unite_rem_cas };
// unite_nd, unite_rem_lock,

/* RemCAS-specific options */
enum SpliceOption {
  split_atomic_one,
  halve_atomic_one,
  splice_simple,
  splice_atomic
};

namespace gbbs {
namespace find_variants {

template<class parent>
inline parent find_naive(parent i, sequence<parent>& parents) {
  while (i != parents[i]) {
    i = parents[i];
  }
  return i;
}

template <class parent>
inline parent find_compress(parent i, sequence<parent>& parents) {
  parent j = i;
  if (parents[j] == j) return j;
  do {
    j = parents[j];
  } while (parents[j] != j);
  parent tmp;
  while ((tmp = parents[i]) > j) {
    parents[i] = j;
    i = tmp;
  }
  return j;
}

template <class parent>
inline parent find_atomic_split(parent i, sequence<parent>& parents) {
  while (1) {
    parent v = parents[i];
    parent w = parents[v];
    if (v == w) {
      return v;
    } else {
      compare_and_swap(&parents[i], v, w);
      // i = its parents
      i = v;
    }
  }
}

template <class parent>
inline parent find_atomic_halve(parent i, sequence<parent>& parents) {
  while (1) {
    parent v = parents[i];
    parent w = parents[v];
    if (v == w) {
      return v;
    } else {
      compare_and_swap(&parents[i], (parent)v, (parent)w);
      // i = its grandparent
      i = parents[i];
    }
  }
}
}  // namespace find_variants

namespace splice_variants {

/* Used in Rem-CAS variants for splice */
template <class parent>
inline parent split_atomic_one(parent i, parent, sequence<parent>& parents) {
  parent v = parents[i];
  parent w = parents[v];
  if (v == w)
    return v;
  else {
    compare_and_swap(&parents[i], v, w);
    i = v;
    return i;
  }
}

/* Used in Rem-CAS variants for splice */
template <class parent>
inline parent halve_atomic_one(parent i, parent, sequence<parent>& parents) {
  parent v = parents[i];
  parent w = parents[v];
  if (v == w)
    return v;
  else {
    compare_and_swap(&parents[i], v, w);
    i = w;
    return i;
  }
}

/* Used in Rem-CAS variants for splice */
template <class parent>
inline parent splice_atomic(parent u, parent v, sequence<parent>& parents) {
  parent z = parents[u];
  compare_and_swap(&parents[u], z, parents[v]);
  return z;
}
}  // namespace splice_variants

namespace unite_variants {

template <class Find, class parent>
struct Unite {
  static constexpr parent UINT_N_MAX = numeric_limits<parent>::max();
  Find& find;
  Unite(Find& find) : find(find) {}

  inline parent operator()(parent u_orig, parent v_orig,
                           sequence<parent>& parents) {
    parent u = u_orig;
    parent v = v_orig;
    while (1) {
      u = find(u, parents);
      v = find(v, parents);
      if (u == v)
        break;
      else if (u > v && parents[u] == u && compare_and_swap(&parents[u], u, v)) {
        return u;
      } else if (v > u && parents[v] == v && compare_and_swap(&parents[v], v, u)) {
        return v;
      }
    }
    return UINT_N_MAX;
  }
};

template <class Splice, class Compress, FindOption find_option, class parent>
struct UniteRemCAS {
  static constexpr parent UINT_N_MAX = numeric_limits<parent>::max();

  Compress& compress;
  Splice& splice;
  UniteRemCAS(Compress& compress, Splice& splice)
      : compress(compress), splice(splice) {}

  inline parent operator()(parent x, parent y, sequence<parent>& parents) {
    parent rx = x;
    parent ry = y;
    while (parents[rx] != parents[ry]) {
      /* link high -> low */
      parent p_ry = parents[ry];
      parent p_rx = parents[rx];
      if (p_rx < p_ry) {
        std::swap(rx, ry);
        std::swap(p_rx, p_ry);
      }
      if (rx == parents[rx] && compare_and_swap(&parents[rx], rx, p_ry)) {
        // if (rx == parents[rx] && (parents[rx] = p_ry)) {
        if constexpr (find_option != find_naive) { /* aka find_none */
          compress(x, parents);
          compress(y, parents);
        }
        return rx;
      } else {
        // failure: locally compress by splicing and try again
        rx = splice(rx, ry, parents);
      }
    }
    return UINT_N_MAX;
  }
};

template <class Find, FindOption find_option, class parent>
struct UniteEarly {
  static constexpr parent UINT_N_MAX = numeric_limits<parent>::max();
  
  Find& find;
  UniteEarly(Find& find) : find(find) {}
  inline parent operator()(parent u, parent v, sequence<parent>& parents) {
    [[maybe_unused]] parent u_orig = u, v_orig = v;
    parent ret = UINT_N_MAX;
    while (u != v) {
      /* link high -> low */
      if (v > u) std::swap(u, v);
      if (parents[u] == u && compare_and_swap(&parents[u], u, v)) {
        ret = u;
        break;
      }
      parent z = parents[u];
      parent w = parents[z];
      compare_and_swap(&parents[u], z, w);
      u = w;
    }
    if constexpr (find_option != find_naive) {
      u = find(u_orig, parents); /* force */
      v = find(v_orig, parents); /* force */
    }
    return ret;
  }
};

}  // namespace unite_variants
}  // namespace gbbs