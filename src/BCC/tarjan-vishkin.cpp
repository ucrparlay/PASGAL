#include "tarjan-vishkin.h"

#include <queue>
#include <type_traits>

#include "graph.h"

constexpr int NUM_ROUND = 5;

template <class Algo, class Graph>
void run(Algo &algo, [[maybe_unused]] const Graph &G,
         [[maybe_unused]] bool verify) {
  using NodeId = typename Graph::NodeId;
  double total_time = 0;
  sequence<NodeId> label;
  for (int i = 0; i <= NUM_ROUND; i++) {
    internal::timer t;
    label = algo.biconnectivity();
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
  algo.get_num_bcc(label);

  ofstream ofs("tarjan-vishkin.tsv", ios_base::app);
  ofs << average_time << '\n';
  ofs.close();

  printf("\n");
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-s] [-v]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-s,\tsymmetrized input graph (required)\n"
            "\t-v,\tverify result\n",
            argv[0]);
    return 0;
  }
  char c;
  char const *input_path = nullptr;
  bool symmetrized = false;
  bool verify = false;
  while ((c = getopt(argc, argv, "i:sv")) != -1) {
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
      default:
        std::cerr << "Error: Unknown option " << optopt << std::endl;
        abort();
    }
  }

  if (!symmetrized) {
    std::cerr << "Error: input graph has to be symmetrized" << std::endl;
    abort();
  }

  printf("Reading graph...\n");
  Graph G;
  G.read_graph(input_path);
  G.symmetrized = symmetrized;

  fprintf(stdout, "Running on %s: |V|=%zu, |E|=%zu, num_round=%d\n", input_path,
          G.n, G.m, NUM_ROUND);

  Tarjan_Vishkin solver(G);
  run(solver, G, verify);

  return 0;
}
