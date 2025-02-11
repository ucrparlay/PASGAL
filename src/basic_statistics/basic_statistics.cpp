#include "parlay/sequence.h"

#include <stack>
#include <algorithm>

#include "graph.h"

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-s] [-v]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-j,\tprint in JSON format\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }
  char c;
  bool return_json = false;
  bool symmetrized false;
  char const *input_path = nullptr;
  while ((c = getopt(argc, argv, "i:j")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      case 'j':
        return_json = true;
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
  // if (!G.symmetrized) {
  G.make_inverse();
  // }

  bool match = true;
  for(size_t i = 0; i < G.n; i++) {
    if(G.offsets[i] != G.in_offsets[i]) {
      match = false;
    }
  }
  for(size_t i = 0; i < G.m; i++) {
    if(G.edges[i] != G.in_edges[i]) {
      match = false;
    }
  }
  if(match) {
    fprintf(stdout, "The graph matched the symmetrized matched.");
  } else {
    fprintf(stdout, "The graph did not matched the symmetrized matched.");
  }
  if(return_json) {
    fprintf(stdout, "{\"vertices_count\": %zu, \"edges_count\": %zu, \"density\": %f, \"avg_degree\": %f}\n",
            G.n, G.m, ((double)G.m * 2) / ((double)G.n * ((double)G.n - 1)), (double)G.m/(double)G.n);
  } else {
    fprintf(stdout, "Running on %s: |V|=%zu, |E|=%zu, density=%f, avg_degree=%f\n",
            input_path, G.n, G.m, ((double)G.m * 2) / ((double)G.n * ((double)G.n - 1)), (double)G.m/(double)G.n);
  }
  return 0;
}
