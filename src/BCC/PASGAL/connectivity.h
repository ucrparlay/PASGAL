#pragma once
#include "ldd.h"
#include "graph.h"
#include "resizable_table.h"
#include "union_find_rules.h"

using namespace std;
using namespace parlay;

template<class NodeId>
struct hash_k {
  uint64_t operator()(const pair<NodeId, NodeId> &k) {
    uint64_t v = k.first;
    v = (v << 32) | (k.second);
    return parlay::hash64(v);
  }
};

template<class NodeId>
NodeId get_max_label(sequence<NodeId> &label) {
  constexpr int NUM_SAMPLES = 1e4;
  size_t n = label.size();
  size_t seed = hash32(n + 1);
  auto samples = sequence<NodeId>::uninitialized(NUM_SAMPLES);
  for (size_t i = 0; i < NUM_SAMPLES; i++) {
    samples[i] = label[hash32(i + seed) % n];
  }
  sort(samples.begin(), samples.end());
  NodeId max_label = 0, max_cnt = 0;
  for (size_t i = 0; i < NUM_SAMPLES;) {
    size_t j = 1;
    while (i + j < NUM_SAMPLES && samples[i + j] == samples[i]) {
      j++;
    }
    if (j > max_cnt) {
      max_cnt = j;
      max_label = samples[i];
    }
    i += j;
  }
  return max_label;
}

template<class Graph, class NodeId = typename Graph::NodeId>
tuple<sequence<NodeId>, sequence<pair<NodeId, NodeId>>> connectivity(
    const Graph &G, double beta,
    function<bool(NodeId, NodeId)> pred = [](NodeId, NodeId) { return true; },
    bool spanning_forest = false) {
  static constexpr NodeId NODE_MAX = numeric_limits<NodeId>::max();
  static constexpr size_t BLOCK_SIZE = 1024;

  LDD ldd_solver(G, pred);
  sequence<NodeId> label;
  sequence<NodeId> parent;
  tie(label, parent) = ldd_solver.ldd(beta, spanning_forest);

  NodeId max_label = get_max_label(label);
  auto find = gbbs::find_variants::find_compress<NodeId>;
  auto splice = gbbs::splice_variants::split_atomic_one<NodeId>;
  auto unite =
      gbbs::unite_variants::UniteRemCAS<decltype(splice), decltype(find),
                                        find_atomic_halve, NodeId>(find, splice);

  auto table = gbbs::resizable_table<pair<NodeId, NodeId>, hash_k<NodeId>>();
  if (spanning_forest) {
    table = gbbs::resizable_table<pair<NodeId, NodeId>, hash_k<NodeId>>(
        G.n, {NODE_MAX, NODE_MAX}, hash_k<NodeId>());
  }
  parallel_for(0, G.n, [&](NodeId i) {
    if (find(label[i], label) != find(max_label, label)) {
      parallel_for(
          G.offsets[i], G.offsets[i + 1],
          [&](size_t j) {
            NodeId v = G.edges[j].v;
            if (pred(i, v)) {
              if (unite(i, v, label) != NODE_MAX) {
                if (spanning_forest) {
                  table.insert({i, v});
                }
              }
            }
          },
          BLOCK_SIZE);
    }
  });

  parallel_for(0, G.n, [&](size_t i) { label[i] = find(label[i], label); });

  sequence<pair<NodeId, NodeId>> tree_edges;
  if (spanning_forest) {
    auto ldd_edges = delayed_seq<pair<NodeId, NodeId>>(G.n, [&](size_t i) {
      return make_pair(i, parent[i]);
    });
    tree_edges = filter(make_slice(ldd_edges), [](const pair<NodeId, NodeId> &a) {
      return a.second != a.first;
    });
    auto union_find_edges = table.entries();
    size_t v1 = tree_edges.size();
    size_t v2 = union_find_edges.size();
    tree_edges.resize(v1 + v2);
    parallel_for(0, v2, [&](size_t i) {
      tree_edges[i + v1] = {get<0>(union_find_edges[i]),
                            get<1>(union_find_edges[i])};
    });
  }
  return {label, tree_edges};
}