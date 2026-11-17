#include "kcore.h"

#include <algorithm>
#include <queue>
#include <vector>

#include "graph.h"
#include "parlay/internal/get_time.h"
#include "parlay/sequence.h"
#include "utils.h"

using namespace std;
using namespace parlay;

constexpr int NUM_ROUND = 5;

template <class Graph, class NodeId = typename Graph::NodeId>
void pal_verifier(const Graph &G, const sequence<NodeId> &act_core) {
  size_t n = G.n;
  size_t num_buckets = 16;
  auto buckets = sequence<hashbag<NodeId>>(num_buckets, hashbag<NodeId>(n));
  auto frontier = sequence<NodeId>::uninitialized(n);
  auto exp_core = sequence<NodeId>::uninitialized(n);
  parallel_for(
      0, n, [&](size_t i) { exp_core[i] = G.offsets[i + 1] - G.offsets[i]; });
  NodeId max_deg = reduce(make_slice(exp_core), maxm<NodeId>());
  NodeId max_core = 0;

  for (NodeId k = 0; k <= max_deg; k += num_buckets) {
    max_deg = 0;
    parallel_for(0, n, [&](size_t i) {
      write_max(&max_deg, exp_core[i]);
      if (exp_core[i] >= k && exp_core[i] < k + num_buckets) {
        buckets[exp_core[i] - k].insert(i);
      }
    });
    if (max_deg < k) {
      break;
    }
    for (NodeId i = 0; i < num_buckets; i++) {
      auto size = buckets[i].pack_into(make_slice(frontier));
      while (size) {
        parallel_for(0, size, [&](size_t j) {
          auto u = frontier[j];
          if (exp_core[u] == k + i) {
            write_max(&max_core, k + i);
            parallel_for(G.offsets[u], G.offsets[u + 1], [&](size_t es) {
              auto v = G.edges[es].v;
              if (exp_core[v] > k + i) {
                auto [id, succeed] =
                    fetch_and_add_bounded(&exp_core[v], -1, k + i);
                id--;
                if (succeed && id + 1 > k + i && id < k + num_buckets) {
                  buckets[id - k].insert(v);
                }
              }
            });
          }
        });
        size = buckets[i].pack_into(make_slice(frontier));
      }
    }
  }
  for (size_t i = 0; i < n; i++) {
    if (exp_core[i] != act_core[i]) {
      printf("exp_core[%zu]: %u while act_core[%zu]: %u\n", i, exp_core[i], i,
             act_core[i]);
    }
    assert(exp_core[i] == act_core[i]);
  }
}

template <class Graph, class NodeId = typename Graph::NodeId>
void verifier(const Graph &G, const sequence<NodeId> &act_core) {
  internal::timer t;
  size_t n = G.n;
  sequence<NodeId> exp_core(n);
  sequence<sequence<NodeId>> buckets(n + 1);
  for (size_t i = 0; i < n; i++) {
    exp_core[i] = G.offsets[i + 1] - G.offsets[i];
    buckets[exp_core[i]].push_back(i);
  }
  NodeId max_core = 0;
  for (NodeId i = 1; i <= n; i++) {
    for (size_t _ = 0; _ < buckets[i].size(); _++) {
      auto u = buckets[i][_];
      if (exp_core[u] == i) {
        max_core = max(max_core, i);
        for (size_t j = G.offsets[u]; j < G.offsets[u + 1]; j++) {
          NodeId v = G.edges[j].v;
          if (exp_core[v] > i) {
            exp_core[v]--;
            buckets[exp_core[v]].push_back(v);
          }
        }
      }
    }
  }
  t.stop();
  printf("Sequential KCore running time: %f\n", t.total_time());
  ofstream ofs("kcore.tsv", ios_base::app);
  ofs << t.total_time() << '\t';
  ofs.close();

  for (size_t i = 0; i < n; i++) {
    if (exp_core[i] != act_core[i]) {
      printf("exp_core[%zu]: %u while act_core[%zu]: %u\n", i, exp_core[i], i,
             act_core[i]);
    }
    assert(exp_core[i] == act_core[i]);
  }
}

template <class Algo, class Graph>
void run(Algo &algo, const Graph &G, bool verify) {
  double total_time = 0;
  using NodeId = typename Graph::NodeId;
  sequence<NodeId> coreness;
  for (int i = 0; i <= NUM_ROUND; i++) {
    internal::timer t;
    coreness = algo.kcore();
    t.stop();
    if (i == 0) {
      printf("Warmup Round: %f\n", t.total_time());
    } else {
      printf("Round %d: %f\n", i, t.total_time());
      total_time += t.total_time();
    }
  }
  double average_time = total_time / NUM_ROUND;
  printf("Average time: %f\n", average_time);
  // printf("Max coreness: %u\n", reduce(coreness, maxm<NodeId>()));

  if (verify) {
    printf("Running verifier...\n");
    internal::timer t;
    auto coreness = algo.kcore();
    t.stop();
    pal_verifier(G, coreness);
  }

  ofstream ofs("kcore.tsv", ios_base::app);
  ofs << average_time << '\n';
  ofs.close();
  printf("\n");
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    fprintf(stderr,
            "Usage: %s [-i input_file] [-s] [-v]\n"
            "Options:\n"
            "\t-i,\tinput file path\n"
            "\t-s,\tsymmetrized input graph\n"
            "\t-v,\tverify result\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }
  char c;
  bool symmetrized = false;
  bool verify = false;
  char const *input_path = nullptr;
  while ((c = getopt(argc, argv, "i:p:a:wsv")) != -1) {
    switch (c) {
      case 'i':
        input_path = optarg;
        break;
      case 's':
        symmetrized = true;
        break;
      case 'v':
        verify = true;
        break;
    }
  }

  printf("Reading graph...\n");
  Graph G;
  G.read_graph(input_path);
  if (!symmetrized) {
    G = make_symmetrized(G);
  }
  G.symmetrized = true;

  printf("Running on %s: |V|=%zu, |E|=%zu, num_round=%d\n", input_path, G.n,
         G.m, NUM_ROUND);

  KCore solver(G);
  run(solver, G, verify);
  return 0;
}
