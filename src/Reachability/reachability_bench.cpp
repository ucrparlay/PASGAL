// Reachability benchmark: sweeps the walk's edge budget (beta), reporting
// time and round count.  Takes many sources per invocation.
#include <getopt.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "graph.h"
#include "reachability.h"

int main(int argc, char *argv[]) {
  char const *input_path = nullptr;
  std::vector<uint32_t> sources;
  bool symmetrized = false;
  bool direction_optimizing = true;
  int measured_rounds = 5;
  std::vector<size_t> betas;

  bool both_directions = false;

  char option;
  while ((option = getopt(argc, argv, "i:sr:n:b:DA")) != -1) {
    switch (option) {
      case 'i': input_path = optarg; break;
      case 's': symmetrized = true; break;
      case 'r': sources.push_back((uint32_t)std::stoul(optarg)); break;
      case 'n': measured_rounds = std::stoi(optarg); break;
      case 'b': betas.push_back((size_t)std::stoul(optarg)); break;
      case 'D': direction_optimizing = false; break;
      case 'A': both_directions = true; break;
      default:
        fprintf(stderr,
                "Usage: %s -i graph [-s] (-r source).. [-n reps] (-b beta).. "
                "[-D | -A]\n  -D disables direction optimization\n"
                "  -A measures both arms, reading the graph once\n",
                argv[0]);
        return EXIT_FAILURE;
    }
  }
  if (input_path == nullptr || sources.empty() || measured_rounds < 1) {
    fprintf(stderr,
            "Usage: %s -i graph [-s] (-r source).. [-n reps] (-b beta).. "
            "[-D | -A]\n",
            argv[0]);
    return EXIT_FAILURE;
  }
  if (betas.empty()) {
    betas = {16, 64, 256, 1024, 4096, 16384, 65536};
  }

  Graph graph;
  graph.read_graph(input_path);
  graph.symmetrized = symmetrized;
  if (!graph.symmetrized) {
    graph.make_inverse();
  }
  for (uint32_t source : sources) {
    if (source >= graph.n) {
      throw std::out_of_range("source vertex is outside the graph");
    }
  }

  // Both arms in one process, so the graph is read and transposed once.
  std::vector<bool> modes;
  if (both_directions) {
    modes = {true, false};
  } else {
    modes = {direction_optimizing};
  }

  Reachability solver(graph);

  for (uint32_t source : sources) {
    printf("REACH_SOURCE\t%u\n", source);
    fflush(stdout);
    // Every beta and both arms must reach the same set; first run is the
    // reference.
    size_t reference = 0;
    bool have_reference = false;
    for (bool diropt : modes) {
      solver.direction_optimizing = diropt;
      for (size_t beta : betas) {
        solver.beta = beta;
        for (int rep = -1; rep < measured_rounds; rep++) {
          parlay::internal::timer timer;
          auto visited = solver.reachability(source);
          timer.stop();
          const size_t reached = parlay::count(visited, (uint8_t)1);
          if (!have_reference) {
            reference = reached;
            have_reference = true;
          } else if (reached != reference) {
            fprintf(stderr,
                    "reachable set changed at beta %zu (%s): %zu vs %zu\n",
                    beta, diropt ? "diropt" : "push", reached, reference);
            return EXIT_FAILURE;
          }
          printf("REACH_TIME\t%s\t%zu\t%s\t%.9f\t%zu\t%zu\t%zu\n",
                 diropt ? "diropt" : "push", beta,
                 rep < 0 ? "warmup" : std::to_string(rep).c_str(),
                 timer.total_time(), solver.rounds(), solver.dense_rounds(),
                 reached);
          fflush(stdout);
        }
      }
    }
  }
  return EXIT_SUCCESS;
}
