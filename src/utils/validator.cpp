#include <queue>
#include <type_traits>

#include "graph.h"

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
  char const *input_path = nullptr;
  bool symmetrized = false;
  while ((c = getopt(argc, argv, "i:ws")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      case 's':
        symmetrized = true;
        break;
      default:
        std::cerr << "Error: Unknown option " << optopt << std::endl;
        abort();
    }
  }

  printf("Reading graph...\n");
  Graph<uint32_t, uint64_t, float> G;
  G.symmetrized = symmetrized;
  G.read_graph(input_path);
  G.validate();
  
  return 0;
}
