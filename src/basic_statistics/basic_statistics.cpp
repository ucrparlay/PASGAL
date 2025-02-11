#include "parlay/sequence.h"

#include <stack>
#include <algorithm>

#include "graph.h"

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-j]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-j,\tprint in JSON format\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }
  char c;
  bool return_json = false;
  char const *input_path = nullptr;
  while ((c = getopt(argc, argv, "i:j")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      case 'j':
        return_json = true;
        break;
    }
  }

  printf("Reading graph...\n");
  Graph G;
  G.read_graph(input_path);

  G.make_inverse();

  size_t max_out_degree = 0, min_out_degree = -1, zero_out_degree_count = 0;
  size_t max_in_degree = 0, min_in_degree = -1, zero_in_degree_count = 0;
  size_t out_degree, in_degree;

  for(size_t i = 0; i < G.n; i++) {
    out_degree = G.offsets[i+1] - G.offsets[i];
    in_degree = G.in_offsets[i+1] - G.in_offsets[i];

    if(out_degree > max_out_degree) {
      max_out_degree = out_degree;
    }
    if(min_out_degree == (size_t)-1 || out_degree < min_out_degree) {
      min_out_degree = out_degree;
    }
    if(out_degree == 0) {
      zero_out_degree_count++;
    }

    if(in_degree > max_in_degree) {
      max_in_degree = in_degree;
    }
    if(min_in_degree == (size_t)-1 || in_degree < min_in_degree) {
      min_in_degree = in_degree;
    }
    if(in_degree == 0) {
      zero_in_degree_count++;
    }
  }

  if(return_json) {
    fprintf(stdout, "{\"vertices_count\": %zu, \"edges_count\": %zu, \"max_out_degree\": %zu, \"min_out_degree\": %zu, \"zero_out_degree_count\": %zu, \"max_in_degree\": %zu, \"min_in_degree\": %zu, \"zero_in_degree_count\": %zu}\n",
            G.n, G.m, max_out_degree, min_out_degree, zero_out_degree_count, max_in_degree, min_in_degree, zero_in_degree_count);
  } else {
    fprintf(stdout, "Running on %s:\n|V| = %zu\n|E| = %zu\nMax Out-degree = %zu\nMin Out-degree = %zu\nZero Out-degree Count = %zu\nMax In-degree = %zu\nMin In-degree = %zu\nZero In-degree Count = %zu\n",
            input_path, G.n, G.m, max_out_degree, min_out_degree, zero_out_degree_count, max_in_degree, min_in_degree, zero_in_degree_count);
  }
  return 0;
}
