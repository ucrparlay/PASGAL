#include "bfs.h"

#include <queue>

#include "graph.h"
#include "seq-bfs.h"

constexpr int NUM_SRC = 5;
constexpr int NUM_ROUND = 5;

template <class Algo, class Graph, class NodeId = typename Graph::NodeId>
void run(Algo &algo, const Graph &G, bool verify, NodeId s) {
  printf("source %-10d\n", s);
  double total_time = 0;
  sequence<NodeId> dist;
  for (int i = 0; i <= NUM_ROUND; i++) {
    internal::timer t;
    dist = algo.bfs(s);
    t.stop();
    if (i == 0) {
      printf("Warmup Round: %f\n", t.total_time());
    } else {
      printf("Round %d: %f\n", i, t.total_time());
      total_time += t.total_time();
    }
  }
  double average_time = total_time / NUM_ROUND;
  printf("Average time: %f\n", average_time);

  ofstream ofs("bfs.tsv", ios_base::app);
  ofs << s << '\t' << average_time << '\n';
  ofs.close();

  if (verify) {
    printf("Running verifier...\n");
    Seq_BFS verifier(G);
    auto exp_dist = verifier.bfs(s);
    assert(dist == exp_dist);
    printf("Passed!\n");
  }
  printf("\n");
}

template <class Algo, class Graph>
void run(Algo &algo, const Graph &G, bool verify) {
  using NodeId = typename Graph::NodeId;
  for (int v = 0; v < NUM_SRC; v++) {
    NodeId s = hash32(v) % G.n;
    run(algo, G, verify, s);
  }
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-s] [-v]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-s,\tsymmetrized input graph\n"
            "\t-v,\tverify result\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }
  char c;
  char const *input_path = nullptr;
  bool symmetrized = false;
  bool verify = false;
  uint32_t source = UINT_MAX;
  while ((c = getopt(argc, argv, "i:svr:")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      case 's':
        symmetrized = true;
        break;
      case 'v':
        verify = true;
        break;
      case 'r':
        source = atol(optarg);
    }
  }

  printf("Reading graph...\n");
  Graph G;
  G.read_graph(input_path);
  G.symmetrized = symmetrized;
  if (!G.symmetrized) {
    G.make_inverse();
  }

  fprintf(stdout, "Running on %s: |V|=%zu, |E|=%zu, num_src=%d, num_round=%d\n",
          input_path, G.n, G.m, NUM_SRC, NUM_ROUND);

  BFS solver(G);
  if (source == UINT_MAX) {
    run(solver, G, verify);
  } else {
    run(solver, G, verify, source);
  }
  return 0;
}
