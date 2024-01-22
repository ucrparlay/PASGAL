#include "dijkstra.h"

#include <queue>
#include <type_traits>

#include "graph.h"

typedef uint32_t NodeId;
typedef uint64_t EdgeId;
#ifdef FLOAT
typedef float EdgeTy;
#else
typedef uint32_t EdgeTy;
#endif

constexpr int NUM_SRC = 10;
constexpr int NUM_ROUND = 1;
constexpr int LOG2_WEIGHT = 18;
constexpr int WEIGHT_RANGE = 1 << LOG2_WEIGHT;

template <class Algo, class Graph, class NodeId = typename Graph::NodeId>
void run(Algo &algo, [[maybe_unused]] const Graph &G, NodeId s) {
  printf("source %-10d\n", s);
  double total_time = 0;
  sequence<EdgeTy> dist;
  for (int i = 0; i <= NUM_ROUND; i++) {
    internal::timer t;
    dist = algo.dijkstra(s);
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

  ofstream ofs("dijkstra.tsv", ios_base::app);
  ofs << s << '\t' << average_time << '\n';
  ofs.close();
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
  Graph<NodeId, EdgeId, EdgeTy> G;
  G.read_graph(input_path);
  G.symmetrized = symmetrized;
  if (!G.weighted) {
    printf("Generating edge weights...\n");
    G.generate_random_weight(1, WEIGHT_RANGE);
  }

  fprintf(stdout, "Running on %s: |V|=%zu, |E|=%zu, num_src=%d, num_round=%d\n",
          input_path, G.n, G.m, NUM_SRC, NUM_ROUND);

  Dijkstra solver(G);
  if (source == UINT_MAX) {
    run(solver, G);
  } else {
    run(solver, G, source);
  }
  return 0;
}
