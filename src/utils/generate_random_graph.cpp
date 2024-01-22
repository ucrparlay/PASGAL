#include "graph.h"

typedef uint32_t NodeId;
typedef uint64_t EdgeId;
typedef uint32_t EdgeTy;

int main(int argc, char* argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-n number of vertices] [-m number of edges] [-o output "
            "file]\n"
            "Options:\n"
            "\t-n,\tnumber of vertices"
            "\t-m,\tnumber of edges"
            "\t-o,\toutput file path\n",
            argv[0]);
    return 0;
  }

  size_t n = 0, m = 0;
  char const* output_path = nullptr;
  char c;
  while ((c = getopt(argc, argv, "n:m:o:")) != -1) {
    switch (c) {
      case 'n':
        n = atoll(optarg);
        break;
      case 'm':
        m = atoll(optarg);
        break;
      case 'o':
        output_path = optarg;
        break;
      default:
        std::cerr << "Error: Unknown option " << optopt << std::endl;
        abort();
    }
  }
  printf("n: %zu, m: %zu\n", n, m);
  Graph<NodeId, EdgeId, EdgeTy> G;
  G.generate_random_graph(n, m);
  G.generate_random_weight(1, 10);
  G = make_symmetrized(G);
  G.write_pbbs_format(output_path);
  return 0;
}
