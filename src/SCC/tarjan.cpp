#include "tarjan.h"

#include <queue>
#include <type_traits>

#include "graph.h"

constexpr int NUM_ROUND = 1;

template <class Algo, class Graph>
void run(Algo &algo, [[maybe_unused]] const Graph &G) {
  using NodeId = typename Graph::NodeId;
  double total_time = 0;
  sequence<NodeId> label;
  for (int i = 0; i <= NUM_ROUND; i++) {
    internal::timer t;
    label = algo.tarjan();
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
  algo.get_num_scc();

  ofstream ofs("tarjan.tsv", ios_base::app);
  ofs << average_time << '\n';
  ofs.close();

  printf("\n");
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file]\n"
            "Options:\n"
            "\t-i,\tinput file path\n",
            argv[0]);
    return 0;
  }
  char c;
  char const *input_path = nullptr;
  while ((c = getopt(argc, argv, "i:")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      default:
        std::cerr << "Error: Unknown option " << optopt << std::endl;
        abort();
    }
  }

  printf("Reading graph...\n");
  Graph G;
  G.read_graph(input_path);
  G.symmetrized = false;

  fprintf(stdout, "Running on %s: |V|=%zu, |E|=%zu, num_round=%d\n", input_path,
          G.n, G.m, NUM_ROUND);

  Tarjan solver(G);
  run(solver, G);

  return 0;
}
