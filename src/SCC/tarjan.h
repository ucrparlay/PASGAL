#pragma once
#include <stack>

#include "graph.h"
using namespace std;
using namespace parlay;

template <class Graph>
class Tarjan {
  using NodeId = typename Graph::NodeId;
  const Graph &G;
  sequence<NodeId> low;
  sequence<NodeId> dfn;
  sequence<NodeId> label;
  stack<NodeId> sk;
  NodeId timestamp;
  size_t num_scc;

  void dfs(NodeId u) {
    low[u] = dfn[u] = ++timestamp;
    sk.push(u);
    for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      NodeId v = G.edges[i].v;
      if (!dfn[v]) {
        dfs(v);
        low[u] = min(low[u], low[v]);
      } else if (!label[v]) {
        low[u] = min(low[u], dfn[v]);
      }
    }
    if (low[u] == dfn[u]) {
      num_scc++;
      while (1) {
        NodeId v = sk.top();
        sk.pop();
        label[v] = num_scc;
        if (v == u) {
          break;
        }
      }
    }
  }

 public:
  Tarjan() = delete;
  Tarjan(const Graph &_G) : G(_G) {
    low = sequence<NodeId>::uninitialized(G.n);
    dfn = sequence<NodeId>::uninitialized(G.n);
    label = sequence<NodeId>::uninitialized(G.n);
  };

  sequence<NodeId> tarjan() {
    for (size_t i = 0; i < G.n; i++) {
      dfn[i] = label[i] = 0;
    }
    sk = stack<NodeId>();
    timestamp = num_scc = 0;
    for (size_t i = 0; i < G.n; i++) {
      if (!dfn[i]) {
        dfs(i);
      }
    }
    return label;
  }

  void get_num_scc() {
    printf("#SCC: %zu\n", num_scc);
    ofstream ofs("tarjan.tsv", ios_base::app);
    ofs << num_scc << '\t';
    ofs.close();
  }
};