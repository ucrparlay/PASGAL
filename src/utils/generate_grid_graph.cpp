#include <climits>
#include <queue>

#include "graph.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "utils.h"

using namespace parlay;
using namespace std;

template <class NodeId = uint32_t, class EdgeId = uint64_t>
Graph<NodeId, EdgeId> generate_grid(size_t row, size_t col, double rate) {
  using Edge = WEdge<NodeId, Empty>;

  size_t n = row * col;
  uint32_t threshold = (uint32_t)(std::numeric_limits<uint32_t>::max() * rate);
  uint32_t threshold2 =
      rate < 0.5 ? threshold * 2 : std::numeric_limits<uint32_t>::max();
  sequence<std::pair<NodeId, NodeId>> edgelist(n * 4);

  constexpr int dx[4] = {1, -1, 0, 0};
  constexpr int dy[4] = {0, 0, 1, -1};
  parallel_for(0, row, [&](size_t x) {
    parallel_for(0, col, [&](size_t y) {
      NodeId u = x * col + y;
      for (int i = 0; i < 4; i++) {
        size_t nx = (x + dx[i] + row) % row, ny = (y + dy[i] + col) % col;
        NodeId v = nx * col + ny;
        assert(u != v);
        uint32_t seed = hash32(u * 4 + i);
        if (seed < threshold) {
          edgelist[u * 4 + i] = make_pair(u, v);
        } else if (seed < threshold2) {
          edgelist[u * 4 + i] = make_pair(v, u);
        } else {
          // A hack to remove the edge since it's a self-loop
          edgelist[u * 4 + i] = make_pair(u, u);
        }
      }
    });
  });
  sort_inplace(make_slice(edgelist));
  auto pred = delayed_seq<bool>(edgelist.size(), [&](size_t i) {
    return edgelist[i].first != edgelist[i].second &&
           (i == 0 || edgelist[i] != edgelist[i - 1]);
  });
  edgelist = pack(edgelist, pred);
  size_t m = edgelist.size();
  Graph<NodeId, EdgeId> G;
  G.n = n;
  G.m = m;
  printf("n: %zu, m: %zu\n", n, m);
  G.offsets = sequence<EdgeId>(n + 1, m);
  G.edges = sequence<Edge>(m);
  parallel_for(0, m, [&](size_t i) {
    G.edges[i].v = edgelist[i].second;
    if (i == 0 || edgelist[i].first != edgelist[i - 1].first) {
      G.offsets[edgelist[i].first] = i;
    }
  });
  scan_inclusive_inplace(make_slice((G.offsets).rbegin(), (G.offsets).rend()),
                         minm<EdgeId>());
  return G;
}

int main(int argc, char* argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-r number of rows] [-c number of columns] [-f sample "
            "rate] [-o output file]\n"
            "Options:\n"
            "\t-r,\tnumber of rows"
            "\t-c,\tnumber of columns"
            "\t-f,\tsample rate"
            "\t-o,\toutput file path\n",
            argv[0]);
    return 0;
  }

  size_t row = 0, col = 0;
  double f = 0;
  char const* output_path = nullptr;
  char c;
  while ((c = getopt(argc, argv, "r:c:f:o:")) != -1) {
    switch (c) {
      case 'r':
        row = atoi(optarg);
        break;
      case 'c':
        col = atoi(optarg);
        break;
      case 'f':
        f = atof(optarg);
      case 'o':
        output_path = optarg;
        break;
      default:
        std::cerr << "Error: Unknown option " << optopt << std::endl;
        abort();
    }
  }
  auto G = generate_grid(row, col, f);
  G.write_binary_format(output_path);
  return 0;
}
