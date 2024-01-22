#include <queue>

#include "BFS/bfs.h"
#include "connectivity.h"
#include "graph.h"

template <class Graph, class NodeId = typename Graph::NodeId>
void get_undirected_diameter(const Graph &G) {
  size_t n = G.n;
  auto cc_ids = get_cc(G);
  sequence<bool> flags(n);
  NodeId diameter = 0;
  parallel_for(0, n, [&](NodeId i) {
    if (compare_and_swap(&flags[cc_ids[i]], false, true)) {
      auto bfs = [&](NodeId s) {
        queue<int> q;
        q.push(s);
        unordered_map<NodeId, NodeId> dist;
        dist[s] = 0;
        NodeId furthest_node = s;
        while (!q.empty()) {
          NodeId u = q.front();
          q.pop();
          for (size_t j = G.offsets[u]; j < G.offsets[u + 1]; j++) {
            NodeId v = G.edges[j].v;
            if (dist.find(v) == dist.end()) {
              dist[v] = dist[u] + 1;
              q.push(v);
              if (dist[v] > dist[furthest_node]) {
                furthest_node = v;
              }
            }
          }
        }
        return make_pair(dist[furthest_node], furthest_node);
      };
      auto furthest_node = bfs(i).second;
      auto dist = bfs(furthest_node).first;
      write_max(&diameter, dist);
    }
  });
  printf("diameter: %u\n", diameter);
  ofstream ofs("diameter.tsv", ios_base::app);
  ofs << diameter << '\n';
  ofs.close();
}

template <class Graph, class NodeId = typename Graph::NodeId>
void get_directed_diameter(const Graph &G) {
  constexpr NodeId DIST_MAX = numeric_limits<NodeId>::max();
  size_t n = G.n;
  NodeId diameter = 0;
  for (size_t i = 0; i < 10000; i++) {
    NodeId t = hash32(i) % n;
    BFS solver(G);
    solver.bfs(t);
    auto dist = solver.bfs(t);
    NodeId furthest_node = t;
    parallel_for(0, n, [&](NodeId j) {
      if (dist[j] != DIST_MAX) {
        write_max(&furthest_node, j,
                  [&](NodeId a, NodeId b) { return dist[a] < dist[b]; });
      }
    });
    write_max(&diameter, dist[furthest_node]);
    dist = solver.bfs(furthest_node);
    parallel_for(0, n, [&](NodeId j) {
      if (dist[j] != DIST_MAX) {
        write_max(&diameter, dist[j]);
      }
    });
  }
  printf("diameter: %u\n", diameter);
  ofstream ofs("diameter.tsv", ios_base::app);
  ofs << diameter << '\n';
  ofs.close();
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-s]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-s,\tsymmetrized input graph\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }
  char c;
  char const *input_path = nullptr;
  bool symmetrized = false;
  while ((c = getopt(argc, argv, "i:s")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      case 's':
        symmetrized = true;
        break;
    }
  }

  printf("Reading graph...\n");
  Graph G;
  G.read_graph(input_path);
  G.symmetrized = symmetrized;
  if (!G.symmetrized) {
    G.make_inverse();
  }

  fprintf(stdout, "Running on %s: |V|=%zu, |E|=%zu\n", input_path, G.n, G.m);

  if (symmetrized) {
    get_undirected_diameter(G);
  } else {
    get_directed_diameter(G);
  }
  return 0;
}
