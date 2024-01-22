#pragma once
#include "connectivity.h"
#include "fast-bcc.h"
#include "resizable_table.h"
#include "spanning_forest.h"
using namespace std;
using namespace parlay;

template <class Graph>
class Tarjan_Vishkin : public BCC<Graph> {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;
  using BCC<Graph>::beta;
  using BCC<Graph>::G;
  using BCC<Graph>::parent;
  using BCC<Graph>::first;
  using BCC<Graph>::last;
  using BCC<Graph>::low;
  using BCC<Graph>::high;
  using BCC<Graph>::euler_tour_tree;
  using BCC<Graph>::tagging;

  static constexpr NodeId NODE_MAX = numeric_limits<NodeId>::max();

 public:
  Tarjan_Vishkin(const Graph &_G) : BCC<Graph>(_G) {}
  sequence<NodeId> biconnectivity() {
    Forest F = spanning_forest(G, beta);
    sequence<NodeId> order = euler_tour_tree(F);
    tagging(F, order);

    sequence<pair<NodeId, NodeId>> edgelist(G.m);
    parallel_for(0, G.n, [&](size_t i) {
      parallel_for(G.offsets[i], G.offsets[i + 1], [&](size_t j) {
        NodeId u = i, v = G.edges[j].v;
        edgelist[j] = make_pair(u, v);
      });
    });
    edgelist = filter(edgelist, [](const pair<NodeId, NodeId> &p) {
      return p.first < p.second;
    });
    edgelist =
        remove_duplicates_ordered(edgelist, less<pair<NodeId, NodeId>>());

    auto parent_edge = sequence<EdgeId>::uninitialized(G.n);
    parallel_for(0, edgelist.size(), [&](size_t i) {
      NodeId u = edgelist[i].first, v = edgelist[i].second;
      if (u < v) {
        if (v == parent[u]) {
          parent_edge[u] = i;
        } else if (u == parent[v]) {
          parent_edge[v] = i;
        }
      }
    });
    // Use NodeId because there is not enough memory to save larger graphs
    auto table = gbbs::resizable_table<pair<NodeId, NodeId>, hash_k<NodeId>>(
        edgelist.size() * 4, {NODE_MAX, NODE_MAX}, hash_k<NodeId>());
    auto check_edge = [&](NodeId u, NodeId v, EdgeId j) {
      if (parent[u] != v && parent[v] != u) {
        if (first[v] < first[u]) {
          table.insert({j, parent_edge[u]});
        }
        if (last[v] <= first[u]) {
          table.insert({parent_edge[v], parent_edge[u]});
        }
      } else {
        if (parent[v] == u) {
          swap(v, u);  // enforce (parent[u] = v)
        }
        if (v != parent[v]) {
          if (low[u] < first[v] || high[u] >= last[v]) {
            table.insert({j, parent_edge[v]});
          }
        }
      }
    };
    parallel_for(0, edgelist.size(), [&](size_t i) {
      NodeId u = edgelist[i].first, v = edgelist[i].second;
      if (u < v) {
        check_edge(u, v, i);
        check_edge(v, u, i);
      }
    });
    auto edges = table.entries();
    auto sym_edges_delayed =
        delayed_seq<pair<NodeId, NodeId>>(edges.size() * 2, [&](size_t i) {
          if (i % 2 == 0) {
            return make_pair(get<0>(edges[i / 2]), get<1>(edges[i / 2]));
          } else {
            return make_pair(get<1>(edges[i / 2]), get<0>(edges[i / 2]));
          }
        });
    auto sym_edges = remove_duplicates_ordered(sym_edges_delayed,
                                               less<pair<NodeId, NodeId>>());
    Graph GA;
    GA.n = edgelist.size();
    GA.m = sym_edges.size();
    using Edge = WEdge<NodeId, Empty>;
    GA.offsets = sequence<EdgeId>(GA.n + 1, GA.m);
    GA.edges = sequence<Edge>::uninitialized(GA.m);
    parallel_for(0, GA.m, [&](size_t i) {
      NodeId u = sym_edges[i].first;
      NodeId v = sym_edges[i].second;
      GA.edges[i].v = v;
      if (i == 0 || sym_edges[i - 1].first != u) {
        GA.offsets[u] = i;
      }
    });
    parlay::scan_inclusive_inplace(
        parlay::make_slice(GA.offsets.rbegin(), GA.offsets.rend()),
        parlay::minm<EdgeId>());
    auto label = get<0>(connectivity(GA, beta));
    return label;
  }

  void get_num_bcc(const sequence<NodeId> &label) {
    size_t num_bcc = remove_duplicates_ordered(label).size();
    printf("#BCC: %zu\n", num_bcc);
    ofstream ofs("tarjan-vishkin.tsv", ios_base::app);
    ofs << num_bcc << '\t';
    ofs.close();
  }
};
