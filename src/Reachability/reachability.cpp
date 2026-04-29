#include "reachability.h"

#include <filesystem>
#include <fstream>
#include <queue>
#include <string>
#include <vector>

#include "BFS/seq-bfs.h"
#include "graph.h"

constexpr int NUM_SRC = 5;
constexpr int NUM_ROUND = 5;

double all_time;

template <class Algo, class Graph, class NodeId = typename Graph::NodeId>
double run(Algo &algo, const Graph &G, bool verify, NodeId s) {
  printf("source %-10d\n", s);
  double total_time = 0;
  sequence<bool> visited;
  for (int i = 0; i <= NUM_ROUND; i++) {
    internal::timer t;
    visited = algo.reachability(s);
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

  all_time += average_time;

  if (verify) {
    constexpr NodeId DIST_MAX = std::numeric_limits<NodeId>::max();
    printf("Running sequential BFS verifier...\n");
    Seq_BFS verifier(G);
    auto exp_dist = verifier.bfs(s);
    for (size_t i = 0; i < G.n; i++) {
      if (visited[i] != (exp_dist[i] != DIST_MAX)) {
        printf("visited[%zu] = %d, exp_dist[%zu] = %d\n", i, visited[i], i,
               exp_dist[i]);
      }
      assert(visited[i] == (exp_dist[i] != DIST_MAX));
    }
    printf("Passed!\n");
  }
  printf("\n");
  return average_time;
}

template <class Algo, class Graph>
void run(Algo &algo, const Graph &G, bool verify,
         const std::string &graph_name, const char *output_path) {
  using NodeId = typename Graph::NodeId;
  std::ofstream ofs(output_path, std::ios_base::app);
  for (int v = 0; v < NUM_SRC; v++) {
    NodeId s = hash32(v) % G.n;
    double average_time = run(algo, G, verify, s);
    ofs << graph_name << '\t' << s << '\t' << average_time << '\t'
        << algo.beta() << '\t' << algo.max_queue_size() << '\t'
        << algo.mode_str() << '\n';
  }
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-o output_tsv] [-s] [-v] [-r source]"
            " [-b beta] [-q max_queue_size] [-t mode]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-o,\tTSV output path\n"
            "\t-s,\tsymmetrized input graph\n"
            "\t-v,\tverify result\n"
            "\t-r,\tsource vertex\n"
            "\t-b,\tbeta (default 2048)\n"
            "\t-q,\tmax_queue_size (default 2000)\n"
            "\t-t,\tthreshold mode: dual | vertex | edge (default vertex)\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }
  char c;
  char const *input_path = nullptr;
  char const *output_path = "reachability.tsv";
  bool symmetrized = false;
  bool verify = false;
  uint32_t source = UINT_MAX;
  size_t beta = 2048;
  size_t max_queue_size = 2000;
  ThresholdMode mode = ThresholdMode::VERTEX;
  while ((c = getopt(argc, argv, "i:o:svr:b:q:t:")) != -1) {
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
      case 'b':
        beta = std::stoul(optarg);
        break;
      case 'q':
        max_queue_size = std::stoul(optarg);
        break;
      case 't': {
        std::string m = optarg;
        if (m == "dual")        mode = ThresholdMode::DUAL;
        else if (m == "vertex") mode = ThresholdMode::VERTEX;
        else if (m == "edge")   mode = ThresholdMode::EDGE;
        else {
          fprintf(stderr, "Unknown threshold mode: %s\n", optarg);
          exit(EXIT_FAILURE);
        }
        break;
      }
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

  Reachability solver(G, beta, max_queue_size, mode);
  fprintf(stdout,
          "Running on %s: |V|=%zu, |E|=%zu, num_src=%d, num_round=%d, "
          "beta=%zu, max_queue_size=%zu, mode=%s\n",
          input_path, G.n, G.m, NUM_SRC, NUM_ROUND, beta, max_queue_size,
          solver.mode_str());
  all_time = 0;
  const std::string graph_name =
      std::filesystem::path(input_path).filename().string();
  if (source == UINT_MAX) {
    run(solver, G, verify, graph_name, output_path);
  } else {
    double average_time = run(solver, G, verify, source);
    std::ofstream ofs(output_path, std::ios_base::app);
    ofs << graph_name << '\t' << source << '\t' << average_time << '\t'
        << solver.beta() << '\t' << solver.max_queue_size() << '\n';
  }
  return 0;
}
