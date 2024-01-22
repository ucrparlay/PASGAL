#pragma once
#include "graph.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "utils.h"

template <class NodeId>
inline NodeId find_compress(NodeId i, parlay::sequence<NodeId> &parents) {
  NodeId j = i;
  if (parents[j] == j) return j;
  do {
    j = parents[j];
  } while (parents[j] != j);
  NodeId tmp;
  while ((tmp = parents[i]) > j) {
    parents[i] = j;
    i = tmp;
  }
  return j;
}

template <class NodeId>
inline bool unite_impl(NodeId u_orig, NodeId v_orig,
                       parlay::sequence<NodeId> &parents) {
  NodeId u = u_orig;
  NodeId v = v_orig;
  while (u != v) {
    u = find_compress(u, parents);
    v = find_compress(v, parents);
    if (u > v && parents[u] == u &&
        atomic_compare_and_swap(&parents[u], u, v)) {
      return true;
    } else if (v > u && parents[v] == v &&
               atomic_compare_and_swap(&parents[v], v, u)) {
      return true;
    }
  }
  return false;
}

// Returns the ids of vertices in the largest connected components.
template <class Graph, class NodeId = typename Graph::NodeId>
parlay::sequence<NodeId> get_cc(const Graph &G) {
  size_t n = G.n;
  parlay::sequence<NodeId> parents(n);
  parlay::parallel_for(0, n, [&](size_t i) { parents[i] = i; });
  parlay::parallel_for(0, n, [&](NodeId u) {
    parlay::parallel_for(
        G.offsets[u], G.offsets[u + 1],
        [&](size_t i) {
          NodeId v = G.edges[i].v;
          unite_impl(u, v, parents);
        },
        1024);
  });
  parlay::parallel_for(
      0, n, [&](NodeId i) { parents[i] = find_compress(i, parents); });
  return parents;
}

// Returns the ids of vertices in the largest connected components.
template <class Graph, class NodeId = typename Graph::NodeId>
parlay::sequence<NodeId> get_largest_cc(const Graph &G) {
  auto parents = get_cc(G);
  auto histogram = parlay::histogram_by_key(parents);
  size_t largest_size = 0;
  parlay::parallel_for(0, histogram.size(), [&](size_t i) {
    write_max(&largest_size, static_cast<size_t>(histogram[i].second));
  });
  NodeId largest_id = 0;
  parlay::parallel_for(0, histogram.size(), [&](size_t i) {
    if (histogram[i].second == largest_size) {
      largest_id = histogram[i].first;
    }
  });
  size_t n = G.n;
  auto pred = parlay::delayed_seq<bool>(
      n, [&](size_t i) { return parents[i] == largest_id; });
  auto largest_cc = parlay::pack_index<NodeId>(pred);
  printf("Size of the largest CC: %zu\n", largest_size);
  return largest_cc;
}
