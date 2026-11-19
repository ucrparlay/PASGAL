#include <queue>

#include "graph.h"
#include "parlay/sequence.h"
using namespace std;
using namespace parlay;

template <class Graph>
class Dijkstra {
  using NodeId = typename Graph::NodeId;
  using EdgeTy = typename Graph::EdgeTy;

  static constexpr EdgeTy DIST_MAX = numeric_limits<EdgeTy>::max();

  const Graph &G;
  sequence<EdgeTy> dist;

 public:
  Dijkstra() = delete;
  Dijkstra(const Graph &_G) : G(_G) {
    dist = sequence<EdgeTy>::uninitialized(G.n);
  }
  sequence<EdgeTy> dijkstra(NodeId s) {
    for (size_t i = 0; i < G.n; i++) {
      dist[i] = DIST_MAX;
    }
    dist[s] = 0;
    priority_queue<pair<EdgeTy, NodeId>, vector<pair<EdgeTy, NodeId>>,
                   greater<pair<EdgeTy, NodeId>>>
        pq;
    pq.push(make_pair(dist[s], s));
    while (!pq.empty()) {
      EdgeTy d;
      NodeId u;
      std::tie(d, u) = pq.top();
      pq.pop();
      if (dist[u] < d) {
        continue;
      }
      for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
        NodeId v = G.edges[i].v;
        EdgeTy w = G.edges[i].w;
        if (dist[v] > dist[u] + w) {
          dist[v] = dist[u] + w;
          pq.push(make_pair(dist[v], v));
        }
      }
    }
    return dist;
  }
};