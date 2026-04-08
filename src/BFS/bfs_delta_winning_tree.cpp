#include "bfs_delta_winning_tree.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <queue>
#include <string>
#include <vector>

#include "graph.h"
#include "seq-bfs.h"

constexpr int NUM_SRC = 5;
constexpr int NUM_ROUND = 5;

template <class Algo, class Graph, class NodeId = typename Graph::NodeId>
double run(Algo &algo, const Graph &G, bool verify, NodeId s) {
  printf("source %-10d\n", s);
  double total_time = 0;
  sequence<NodeId> dist;
  for (int i = 0; i < NUM_ROUND; i++) {
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
  constexpr int measured_rounds = NUM_ROUND - 1;
  const double average_time = total_time / measured_rounds;
  printf("Average time: %f\n", average_time);

  if (verify) {
    printf("Running verifier...\n");
    Seq_BFS verifier(G);
    auto exp_dist = verifier.bfs(s);
    if (dist != exp_dist) {
      int cnt = 0;
      for (size_t i = 0; i < G.n; i++) {
        if (dist[i] != exp_dist[i]) {
          printf("dist[%zu] = %d, exp_dist[%zu] = %d\n", i, dist[i], i,
                 exp_dist[i]);
          fflush(stdout);
          if (cnt++ > 10) {
            break;
          }
        }
      }
    }
    assert(dist == exp_dist);
    printf("Passed!\n");
  }
  printf("\n");
  return average_time;
}

template <class Algo, class Graph>
double run(Algo &algo, const Graph &G, bool verify) {
  using NodeId = typename Graph::NodeId;
  double log_sum = 0.0;
  for (int v = 0; v < NUM_SRC; v++) {
    NodeId s = hash32(v) % G.n;
    double average_time = run(algo, G, verify, s);
    log_sum += std::log(average_time);
  }
  double geomean = std::exp(log_sum / NUM_SRC);
  printf("Geomean: %f\n\n", geomean);
  return geomean;
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-o output_tsv] [-s] [-v] [-r source]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-o,\tTSV output path\n"
            "\t-s,\tsymmetrized input graph\n"
            "\t-v,\tverify result\n"
            "\t-r,\tsource vertex\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }
  char c;
  char const *input_path = nullptr;
  char const *output_path = "bfs_delta_winning_tree.tsv";
  bool symmetrized = false;
  bool verify = false;
  uint32_t source = UINT_MAX;
  while ((c = getopt(argc, argv, "i:o:svr:")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      case 'o':
        output_path = optarg;
        break;
      case 's':
        symmetrized = true;
        break;
      case 'v':
        verify = true;
        break;
      case 'r':
        source = atol(optarg);
        break;
      default:
        fprintf(stderr, "Unknown option: -%c\n", c);
        exit(EXIT_FAILURE);
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

  std::vector<double> step_times;
  step_times.reserve(11);
  for (int step = 1; step <= (1 << 10); step *= 2) {
    printf("step: %d\n", step);
    BFS solver(G, step);
    if (source == UINT_MAX) {
      step_times.push_back(run(solver, G, verify));
    } else {
      step_times.push_back(run(solver, G, verify, source));
    }
  }

  const std::string graph_name = std::filesystem::path(input_path).filename().string();
  std::ofstream ofs(output_path, std::ios_base::app);
  ofs << graph_name;
  for (double time : step_times) {
    ofs << char(9) << time;
  }
  ofs << char(10);
  return 0;
}
