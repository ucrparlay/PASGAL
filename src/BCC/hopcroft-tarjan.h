#pragma once
#include <stack>

#include "graph.h"
using namespace std;
using namespace parlay;

template <class Graph>
class Hopcroft_Tarjan {
 protected:
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr NodeId NODE_MAX = numeric_limits<NodeId>::max();

  const Graph &G;
  sequence<NodeId> low;
  sequence<NodeId> dfn;
  sequence<NodeId> label;
  stack<pair<NodeId, NodeId>> sk;
  NodeId idx;
  size_t num_bcc;

  void dfs(NodeId u, NodeId f = NODE_MAX) {
    dfn[u] = low[u] = ++idx;
    for (size_t j = G.offsets[u]; j < G.offsets[u + 1]; j++) {
      auto v = G.edges[j].v;
      if (!dfn[v]) {
        sk.push({u, v});
        dfs(v, u);
        low[u] = min(low[u], low[v]);
        if (low[v] >= dfn[u]) {
          num_bcc++;
          while (true) {
            auto [f, s] = sk.top();
            sk.pop();
            label[f] = label[s] = num_bcc;
            if (f == u && s == v) {
              break;
            }
          }
        }
      } else if (v != f) {
        low[u] = min(low[u], dfn[v]);
      }
    }
  }

 public:
  Hopcroft_Tarjan() = delete;
  Hopcroft_Tarjan(const Graph &_G) : G(_G) {
    low = sequence<NodeId>::uninitialized(G.n);
    dfn = sequence<NodeId>::uninitialized(G.n);
    label = sequence<NodeId>::uninitialized(G.n);
  }

  sequence<NodeId> biconnectivity() {
    idx = num_bcc = 0;
    for (size_t i = 0; i < G.n; i++) {
      dfn[i] = 0;
    }
    while(!sk.empty()) {
      sk.pop();
    }
    for (size_t i = 0; i < G.n; i++) {
      if (!dfn[i]) {
        dfs(i);
      }
    }
    return label;
  }

  void get_num_bcc() {
    printf("#BCC: %zu\n", num_bcc);
    ofstream ofs("hopcroft-tarjan.tsv", ios_base::app);
    ofs << num_bcc << '\t';
    ofs.close();
  }
};