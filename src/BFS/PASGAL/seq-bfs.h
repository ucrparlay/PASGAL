#include <queue>

#include "graph.h"
#include "parlay/sequence.h"
using namespace std;
using namespace parlay;

template <class Graph>
class Seq_BFS {
  using NodeId = typename Graph::NodeId;

  static constexpr NodeId DIST_MAX = numeric_limits<NodeId>::max();

  const Graph &G;
  sequence<NodeId> dist;

 public:
  Seq_BFS() = delete;
  Seq_BFS(const Graph &_G) : G(_G) {
    dist = sequence<NodeId>::uninitialized(G.n);
  }
  sequence<NodeId> bfs(NodeId s) {
    for (size_t i = 0; i < G.n; i++) {
      dist[i] = DIST_MAX;
    }
    dist[s] = 0;
    queue<int> q;
    q.push(s);
    while (!q.empty()) {
      NodeId u = q.front();
      q.pop();
      for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
        NodeId v = G.edges[i].v;
        if (dist[v] == DIST_MAX) {
          dist[v] = dist[u] + 1;
          q.push(v);
        }
      }
    }
    return dist;
  }
};
