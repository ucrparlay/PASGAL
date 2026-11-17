#include "bfs.h"

#include <queue>

#include "connectivity.h"
#include "graph.h"
#include "seq-bfs.h"

constexpr int NUM_SRC = 10;
constexpr int NUM_ROUND = 5;

template <class Func>
void run(Func &&f) {
  double total_time = 0;
  for (int i = 0; i <= NUM_ROUND; i++) {
    internal::timer tm;
    f();
    tm.stop();
    if (i == 0) {
      printf("Warmup Round: %f\n", tm.total_time());
    } else {
      printf("Round %d: %f\n", i, tm.total_time());
      total_time += tm.total_time();
    }
  }
  double average_time = total_time / NUM_ROUND;
  printf("Average time: %f\n", average_time);

  ofstream ofs("bfs.tsv", ios_base::app);
  ofs << average_time << '\t';
  ofs.close();
}

template <class Graph, class NodeId = typename Graph::NodeId>
void run(const Graph &G, const sequence<NodeId> &largest_cc) {
  size_t m = largest_cc.size();
  for (int v = 0; v < NUM_SRC; v++) {
    NodeId s = largest_cc[hash32(v) % m];
    sequence<NodeId> dist1, dist2;

    {
      BFS solver(G);
      run([&]() { dist1 = solver.bfs(s); });
    }

    {
      Seq_BFS solver(G);
      run([&]() { dist2 = solver.bfs(s); });
    }

    assert(dist1 == dist2);
    ofstream ofs("bfs.tsv", ios_base::app);
    ofs << s << '\n';
    printf("\n");
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
  bool symmetrized = false;
  char const *input_path = nullptr;
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

  fprintf(stdout, "Running on %s: |V|=%zu, |E|=%zu, num_src=%d, num_round=%d\n",
          input_path, G.n, G.m, NUM_SRC, NUM_ROUND);

  auto largest_cc = get_largest_cc(G);
  run(G, largest_cc);
  return 0;
}
