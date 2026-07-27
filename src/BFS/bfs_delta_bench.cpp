#include "bfs_delta.h"

#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "graph.h"

template <class Sequence>
void write_distances(const std::string &path, uint64_t num_vertices,
                     uint32_t source, const Sequence &distances) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    throw std::runtime_error("could not open distance output: " + path);
  }

  const char magic[8] = {'B', 'F', 'S', 'D', 'S', 'T', '0', '1'};
  uint64_t reachable = 0;
  const uint32_t reserved = 0;
  out.write(magic, sizeof(magic));
  out.write(reinterpret_cast<const char *>(&num_vertices),
            sizeof(num_vertices));
  const std::streampos reachable_pos = out.tellp();
  out.write(reinterpret_cast<const char *>(&reachable), sizeof(reachable));
  out.write(reinterpret_cast<const char *>(&source), sizeof(source));
  out.write(reinterpret_cast<const char *>(&reserved), sizeof(reserved));

  for (uint64_t vertex = 0; vertex < num_vertices; vertex++) {
    const uint32_t distance = distances[vertex];
    if (distance != std::numeric_limits<uint32_t>::max()) {
      const uint32_t vertex_id = static_cast<uint32_t>(vertex);
      out.write(reinterpret_cast<const char *>(&vertex_id), sizeof(vertex_id));
      out.write(reinterpret_cast<const char *>(&distance), sizeof(distance));
      reachable++;
    }
  }
  out.seekp(reachable_pos);
  out.write(reinterpret_cast<const char *>(&reachable), sizeof(reachable));
  if (!out) {
    throw std::runtime_error("failed while writing distance output: " + path);
  }
}

template <class Sequence>
void verify_distances(const Sequence &actual, const Sequence &expected,
                      uint32_t delta, int repetition) {
  if (actual.size() != expected.size()) {
    throw std::runtime_error("bfs_delta distance array size changed");
  }
  for (size_t vertex = 0; vertex < actual.size(); vertex++) {
    if (actual[vertex] != expected[vertex]) {
      throw std::runtime_error(
          "bfs_delta distances differ at delta " + std::to_string(delta) +
          ", repetition " + std::to_string(repetition) + ", vertex " +
          std::to_string(vertex) + ": expected " +
          std::to_string(expected[vertex]) + ", got " +
          std::to_string(actual[vertex]));
    }
  }
}

int main(int argc, char *argv[]) {
  char const *input_path = nullptr;
  // Repeatable: several sources per invocation amortize the graph read and,
  // for a directed graph, the transpose.  -o is given once per -r, in the
  // same order.
  std::vector<uint32_t> sources;
  std::vector<std::string> distance_paths;
  bool symmetrized = false;
  int measured_rounds = 10;

  char option;
  while ((option = getopt(argc, argv, "i:o:sr:n:")) != -1) {
    switch (option) {
      case 'i':
        input_path = optarg;
        break;
      case 'o':
        distance_paths.emplace_back(optarg);
        break;
      case 's':
        symmetrized = true;
        break;
      case 'r':
        sources.push_back(static_cast<uint32_t>(std::stoul(optarg)));
        break;
      case 'n':
        measured_rounds = std::stoi(optarg);
        break;
      default:
        fprintf(stderr,
                "Usage: %s -i graph (-r source -o distances).. [-s] "
                "[-n rounds]\n",
                argv[0]);
        return EXIT_FAILURE;
    }
  }
  if (input_path == nullptr || sources.empty() || measured_rounds < 1 ||
      distance_paths.size() != sources.size()) {
    fprintf(stderr,
            "Usage: %s -i graph (-r source -o distances).. [-s] "
            "[-n rounds]\n  -r and -o repeat together, one -o per -r\n",
            argv[0]);
    return EXIT_FAILURE;
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

  for (size_t index = 0; index < sources.size(); index++) {
    const uint32_t source = sources[index];
    // Marks which source the timings that follow belong to.  Readers that
    // pass a single source can ignore it.
    printf("BFS_BENCH_SOURCE\t%u\n", source);
    fflush(stdout);
    parlay::sequence<uint32_t> reference_distances;
    for (uint32_t delta = 1; delta <= (1U << 10); delta *= 2) {
      BFS solver(graph, delta);
      for (int repetition = -1; repetition < measured_rounds; repetition++) {
        solver.prepare(source);
        internal::timer timer;
        solver.run_prepared();
        timer.stop();
        printf("BFS_BENCH_TIME\t%u\t%s\t%.9f\n", delta,
               repetition < 0 ? "warmup"
                              : std::to_string(repetition).c_str(),
               timer.total_time());
        fflush(stdout);
        if (delta == 1 && repetition < 0) {
          reference_distances = solver.distances();
          write_distances(distance_paths[index], graph.n, source,
                          reference_distances);
        } else {
          verify_distances(solver.distances(), reference_distances, delta,
                           repetition);
        }
      }
    }
  }
  return EXIT_SUCCESS;
}
