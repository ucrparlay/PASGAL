#include "seq-bfs.h"

#include <queue>
#include <type_traits>

#include "graph.h"

constexpr int NUM_SRC = 5;
constexpr int NUM_ROUND = 1;

template <class Algo, class Graph, class NodeId = typename Graph::NodeId>
void run(Algo &algo, [[maybe_unused]] const Graph &G, NodeId s) {
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

  ofstream ofs("seq-bfs.tsv", ios_base::app);
  ofs << s << '\t' << average_time << '\n';
  ofs.close();

  printf("\n");
}

template <class Algo, class Graph>
void run(Algo &algo, const Graph &G) {
  using NodeId = typename Graph::NodeId;
  for (int v = 0; v < NUM_SRC; v++) {
    NodeId s = hash32(v) % G.n;
    run(algo, G, s);
  }
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-s]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-s,\tsymmetrized input graph\n",
            argv[0]);
    return 0;
  }
  char c;
  char const *input_path = nullptr;
  bool symmetrized = false;
  uint32_t source = UINT_MAX;
  while ((c = getopt(argc, argv, "i:sr:")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      case 's':
        symmetrized = true;
        break;
      case 'r':
        source = atol(optarg);
        break;
      default:
        std::cerr << "Error: Unknown option " << optopt << std::endl;
        abort();
    }
  }

  printf("Reading graph...\n");
  Graph G;
  G.read_graph(input_path);
  G.symmetrized = symmetrized;

  fprintf(stdout, "Running on %s: |V|=%zu, |E|=%zu, num_src=%d, num_round=%d\n",
          input_path, G.n, G.m, NUM_SRC, NUM_ROUND);

  Seq_BFS solver(G);
  if (source == UINT_MAX) {
    run(solver, G);
  } else {
    run(solver, G, source);
  }
  return 0;
}
