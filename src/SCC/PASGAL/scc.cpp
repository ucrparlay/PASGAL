#include "scc.h"
#include "reach.h"  // remove later, included in scc.h
#include "multi_reach.h"
#include "parlay/sequence.h"


#include <stack>
#include <algorithm>

#include "graph.h"

constexpr int NUM_ROUND = 10;


int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-s] [-v]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-v,\tverify result\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }
  char c;
  // bool verify = false;
  char const *input_path = nullptr;
  while ((c = getopt(argc, argv, "i:p:a:wsv")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      // case 'v':
      //   verify = true;
        // break;
    }
  }

  printf("Reading graph...\n");
  Graph G;
  // Graph G;
  // G.read_graph(input_path);
  G.read_graph(input_path);
  auto GT = Transpose(G);
  fprintf(stdout, "Running on %s: |V|=%zu, |E|=%zu, num_round=%d\n",
          input_path, G.n, G.m, NUM_ROUND);
  
  SCC scc_solver(G, GT);
  parlay::sequence<uint64_t> scc_labels=parlay::sequence<uint64_t>(G.n);
  parlay::internal::timer t("SCC");
  scc_solver.scc(scc_labels);
  int r = 5;
  t.start();
  for (int i = 0; i<r;i++){
    scc_solver.scc(scc_labels);
    t.next("scc");
  }
  printf("average scc cost: %f\n", t.total_time()/r);
  scc_solver.status(scc_labels);
  return 0;
}