#include "graph.h"

typedef uint32_t NodeId;
typedef uint64_t EdgeId;
typedef float EdgeTy;

int main(int argc, char* argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-o output file]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-o,\toutput file path\n",
            argv[0]);
    return 0;
  }

  char const* input_path = nullptr;
  char const* output_path = nullptr;
  char c;
  while ((c = getopt(argc, argv, "i:o:")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      case 'o':
        output_path = optarg;
        break;
      default:
        std::cerr << "Error: Unknown option " << optopt << std::endl;
        abort();
    }
  }
  printf("Reading graph...\n");
  Graph<NodeId, EdgeId, EdgeTy> G;
  G.read_graph(input_path);
  G = make_symmetrized(G);
  //G.write_pbbs_format(output_path);
  G.write_binary_format(output_path);
  return 0;
}
