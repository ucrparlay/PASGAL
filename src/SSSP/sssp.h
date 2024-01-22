#pragma once
#include <climits>

#include "graph.h"
#include "hashbag.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"
#include "utils.h"

using namespace std;
using namespace parlay;

enum Algorithm { rho_stepping = 0, delta_stepping, bellman_ford };

template <class Graph>
class SSSP {
 protected:
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;
  using EdgeTy = typename Graph::EdgeTy;

  static constexpr EdgeTy DIST_MAX = numeric_limits<EdgeTy>::max();
  static constexpr size_t LOCAL_QUEUE_SIZE = 128;
  static constexpr size_t BLOCK_SIZE = 1024;
  static constexpr size_t NUM_SAMPLES = 1024;
  static constexpr size_t SPARSE_TH = 100;
  static constexpr size_t GROWTH_FACTOR = 10;

  const Graph &G;
  bool sparse;
  size_t frontier_size;
  hashbag<NodeId> bag;
  sequence<EdgeTy> dist;
  sequence<NodeId> frontier;
  sequence<atomic<bool>> in_frontier;
  sequence<atomic<bool>> in_next_frontier;

  virtual void init() = 0;
  virtual EdgeTy get_threshold() = 0;

  void add_to_frontier(NodeId v) {
    if (sparse) {
      if (!in_frontier[v] &&
          compare_and_swap(&in_next_frontier[v], false, true)) {
        bag.insert(v);
      }
    } else {  // dense
      if (!in_frontier[v] && !in_next_frontier[v]) {
        in_next_frontier[v] = true;
      }
    }
  }

  size_t estimate_size() {
    static uint32_t seed = G.n;
    size_t hits = 0;
    for (size_t i = 0; i < NUM_SAMPLES; i++) {
      NodeId u = hash32(seed) % G.n;
      if (in_frontier[u]) {
        hits++;
      }
      seed++;
    }
    return hits * G.n / NUM_SAMPLES;
  }

  inline void visit_neighbors_sequential(NodeId u, NodeId *local_queue,
                                         size_t &rear) {
    if (G.symmetrized) {
      EdgeTy min_dist = dist[u];
      for (EdgeId i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
        NodeId v = G.edges[i].v;
        EdgeTy w = G.edges[i].w;
        if (dist[v] != DIST_MAX) {
          min_dist = min(min_dist, dist[v] + w);
        }
      }
      write_min(&dist[u], min_dist);
    }
    for (EdgeId i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      NodeId v = G.edges[i].v;
      EdgeTy w = G.edges[i].w;
      if (write_min(&dist[v], dist[u] + w)) {
        if (rear < LOCAL_QUEUE_SIZE) {
          local_queue[rear++] = v;
        } else {
          add_to_frontier(v);
        }
      }
    }
  }

  inline void visit_neighbors_parallel(NodeId u) {
    blocked_for(G.offsets[u], G.offsets[u + 1], BLOCK_SIZE,
                [&](size_t, size_t start, size_t end) {
                  if (G.symmetrized) {
                    EdgeTy min_dist = dist[u];
                    for (EdgeId i = start; i < end; i++) {
                      NodeId v = G.edges[i].v;
                      EdgeTy w = G.edges[i].w;
                      if (dist[v] != DIST_MAX) {
                        min_dist = min(min_dist, dist[v] + w);
                      }
                    }
                    if (write_min(&dist[u], min_dist)) {
                      add_to_frontier(u);
                    }
                  }
                  for (EdgeId i = start; i < end; i++) {
                    NodeId v = G.edges[i].v;
                    EdgeTy w = G.edges[i].w;
                    if (write_min(&dist[v], dist[u] + w)) {
                      add_to_frontier(v);
                    }
                  }
                });
  }

  size_t sparse_relax() {
    constexpr bool use_local_queue = true;

    EdgeTy threshold = get_threshold();
    parallel_for(
        0, frontier_size,
        [&](size_t i) {
          NodeId f = frontier[i];
          in_frontier[f] = false;
          if (dist[f] > threshold) {
            add_to_frontier(f);
          } else {
            if (use_local_queue) {
              NodeId local_queue[LOCAL_QUEUE_SIZE];
              size_t front = 0, rear = 0;
              local_queue[rear++] = f;
              while (front < rear) {
                NodeId u = local_queue[front++];
                if (dist[u] > threshold) {
                  add_to_frontier(u);
                  continue;
                }
                size_t deg = G.offsets[u + 1] - G.offsets[u];
                if (deg < 10 * BLOCK_SIZE) {
                  visit_neighbors_sequential(u, local_queue, rear);
                } else {
                  visit_neighbors_parallel(u);
                }
              }
            } else {
              visit_neighbors_parallel(f);
            }
          }
        },
        1);
    swap(in_frontier, in_next_frontier);
    return bag.pack_into(make_slice(frontier));
  }

  size_t dense_relax() {
    while (estimate_size() >= G.n / SPARSE_TH) {
      EdgeTy threshold = get_threshold();
      parallel_for(
          0, G.n,
          [&](NodeId u) {
            if (in_frontier[u]) {
              in_frontier[u] = false;
              if (dist[u] > threshold) {
                add_to_frontier(u);
              } else {
                visit_neighbors_parallel(u);
              }
            }
          },
          1);
      swap(in_frontier, in_next_frontier);
    }
    return count(in_frontier, true);
  }

  void sparse2dense() {}

  void dense2sparse() {
    auto identity = delayed_seq<NodeId>(G.n, [&](NodeId i) { return i; });
    pack_into_uninitialized(identity, in_frontier, frontier);
  }

 public:
  SSSP() = delete;
  SSSP(const Graph &_G) : G(_G), bag(G.n) {
    dist = sequence<EdgeTy>::uninitialized(G.n);
    frontier = sequence<NodeId>::uninitialized(G.n);
    in_frontier = sequence<atomic<bool>>::uninitialized(G.n);
    in_next_frontier = sequence<atomic<bool>>::uninitialized(G.n);
  }

  sequence<EdgeTy> sssp(NodeId s) {
    if (!G.weighted) {
      fprintf(stderr, "Error: Input graph is unweighted\n");
      exit(EXIT_FAILURE);
    }

    init();
    parallel_for(0, G.n, [&](NodeId i) {
      dist[i] = DIST_MAX;
      in_frontier[i] = in_next_frontier[i] = false;
    });
    assert(bag.pack_into(make_slice(frontier)) == 0);

    frontier_size = 1;
    dist[s] = 0;
    frontier[0] = s;
    in_frontier[s] = true;
    sparse = true;

    // int round = 0;
    while (frontier_size) {
      // printf("Round %d: %s, size: %zu, ", round++, sparse ? "sparse" :
      // "dense", frontier_size);
      // internal::timer t;
      if (sparse) {
        frontier_size = sparse_relax();
      } else {
        frontier_size = dense_relax();
      }
      // printf("relax: %f, ", t.next_time());
      bool next_sparse = (frontier_size < G.n / SPARSE_TH) ? true : false;
      if (sparse && !next_sparse) {
        sparse2dense();
      } else if (!sparse && next_sparse) {
        dense2sparse();
      }
      // printf("pack: %f\n", t.next_time());
      sparse = next_sparse;
    }
    return dist;
  }
};

template <class Graph>
class Rho_Stepping : public SSSP<Graph> {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;
  using EdgeTy = typename Graph::EdgeTy;
  using SSSP<Graph>::frontier_size;
  using SSSP<Graph>::sparse;
  using SSSP<Graph>::dist;
  using SSSP<Graph>::frontier;
  using SSSP<Graph>::G;
  using SSSP<Graph>::in_frontier;

  static constexpr EdgeTy DIST_MAX = numeric_limits<EdgeTy>::max();
  static constexpr size_t NUM_SAMPLES = 1024;

  size_t rho;
  uint32_t seed;
  void init() override { seed = 0; }
  EdgeTy get_threshold() override {
    if (frontier_size <= rho) {
      if (sparse) {
        auto _dist = delayed_seq<EdgeTy>(
            frontier_size, [&](size_t i) { return dist[frontier[i]]; });
        return *max_element(_dist);
      } else {
        return DIST_MAX;
      }
    }
    EdgeTy sample_dist[NUM_SAMPLES + 1];
    for (size_t i = 0; i <= NUM_SAMPLES; i++) {
      if (sparse) {
        NodeId v = frontier[hash32(seed + i) % frontier_size];
        sample_dist[i] = dist[v];
      } else {
        NodeId v = hash32(seed + i) % G.n;
        if (in_frontier[v]) {
          sample_dist[i] = dist[v];
        } else {
          sample_dist[i] = DIST_MAX;
        }
      }
    }
    seed += NUM_SAMPLES + 1;
    size_t id = 1.0 * rho / frontier_size * NUM_SAMPLES;
    sort(sample_dist, sample_dist + NUM_SAMPLES + 1);
    return sample_dist[id];
  }

 public:
  Rho_Stepping(const Graph &_G, size_t _rho = 1 << 20)
      : SSSP<Graph>(_G), rho(_rho) {}
};

template <class Graph>
class Delta_Stepping : public SSSP<Graph> {
  using EdgeTy = typename Graph::EdgeTy;

  EdgeTy delta;
  EdgeTy thres;

  void init() override { thres = 0; }
  EdgeTy get_threshold() override {
    thres += delta;
    return thres;
  }

 public:
  Delta_Stepping(const Graph &_G, EdgeTy _delta = 1 << 15)
      : SSSP<Graph>(_G), delta(_delta) {}
};

template <class Graph>
class Bellman_Ford : public SSSP<Graph> {
  using NodeId = typename Graph::NodeId;
  using EdgeTy = typename Graph::EdgeTy;

  static constexpr EdgeTy DIST_MAX = numeric_limits<EdgeTy>::max();

  void init() override {}
  EdgeTy get_threshold() override { return DIST_MAX; }

 public:
  Bellman_Ford(const Graph &_G) : SSSP<Graph>(_G) {}
};
