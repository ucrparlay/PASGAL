#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>

#include "graph.h"
#include "hashbag.h"
#include "parlay/parallel.h"
#include "parlay/random.h"
#include "parlay/sequence.h"
#include "parlay/primitives.h"

using namespace std;
using namespace parlay;

template <class Graph>
class LDD {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr NodeId NODE_MAX = numeric_limits<NodeId>::max();
  static constexpr size_t LOCAL_QUEUE_SIZE = 1024;
  static constexpr size_t BLOCK_SIZE = 1024;
  static constexpr size_t NUM_SAMPLES = 1024;

 private:
  const Graph& G;
  sequence<NodeId> frontier;
  sequence<bool> in_frontier;
  sequence<bool> in_next_frontier;
  hashbag<NodeId> bag;
  bool sparse;
  size_t frontier_size;
  size_t threshold;
  function<bool(NodeId, NodeId)> pred;

  size_t sparse_update(sequence<NodeId>& label, sequence<NodeId>& parent) {
    parallel_for(
        0, frontier_size,
        [&](size_t i) {
          NodeId f = frontier[i];
          size_t deg_f = G.offsets[f + 1] - G.offsets[f];
          if (deg_f > BLOCK_SIZE) {
            parallel_for(
                G.offsets[f], G.offsets[f + 1],
                [&](size_t j) {
                  NodeId v = G.edges[j].v;
                  if (pred(f, v)) {
                    if (compare_and_swap(&label[v], NODE_MAX, label[f])) {
                      if (parent.size()) {
                        parent[v] = f;
                      }
                      bag.insert(v);
                    }
                  }
                },
                BLOCK_SIZE);
          } else {
            NodeId local_queue[LOCAL_QUEUE_SIZE];
            size_t head = 0, tail = 0;
            local_queue[tail++] = f;
            while (head < tail && tail != LOCAL_QUEUE_SIZE) {
              NodeId u = local_queue[head++];
              size_t deg_u = G.offsets[u + 1] - G.offsets[u];
              if (deg_u > BLOCK_SIZE) {
                bag.insert(u);
              } else {
                for (size_t j = G.offsets[u]; j < G.offsets[u + 1]; j++) {
                  NodeId v = G.edges[j].v;
                  if (pred(u, v)) {
                    if (compare_and_swap(&label[v], NODE_MAX, label[u])) {
                      if (parent.size()) {
                        parent[v] = u;
                      }
                      if (tail < LOCAL_QUEUE_SIZE) {
                        local_queue[tail++] = v;
                      } else {
                        bag.insert(v);
                      }
                    }
                  }
                }
              }
            }
            for (size_t j = head; j < tail; j++) {
              bag.insert(local_queue[j]);
            }
          }
        },
        1);
    return bag.pack_into(make_slice(frontier));
  }

  size_t dense_update(sequence<NodeId>& label, sequence<NodeId>& parent) {
    parallel_for(
        0, G.n,
        [&](NodeId i) {
          in_next_frontier[i] = false;
          if (label[i] == NODE_MAX) {
            const auto& offsets = G.symmetrized ? G.offsets : G.in_offsets;
            const auto& edges = G.symmetrized ? G.edges : G.in_edges;
            for (size_t j = offsets[i]; j < offsets[i + 1]; j++) {
              NodeId v = edges[j].v;
              if (pred(i, v)) {
                if (in_frontier[v]) {
                  if (parent.size()) {
                    parent[i] = v;
                  }
                  label[i] = label[v];
                  in_frontier[i] = in_next_frontier[i] = true;
                  break;
                }
              }
            }
          }
        },
        BLOCK_SIZE);

    swap(in_frontier, in_next_frontier);
    auto ret = count(in_frontier, true);
    return ret;
  }

  EdgeId dense_sample(NodeId seed) {
    constexpr int NUM_SAMPLES = 50;
    NodeId count = 0;
    EdgeId out_edges = 0;
    int i = 0;
    while (count < NUM_SAMPLES) {
      i++;
      NodeId index = hash32(seed + i) % G.n;
      if (in_frontier[index]) {
        count++;
        out_edges += G.offsets[index + 1] - G.offsets[index];
      }
    }
    return frontier_size * (out_edges / count);
  }

  void sparse2dense() {
    parallel_for(0, G.n, [&](size_t i) { in_frontier[i] = false; });
    parallel_for(0, frontier_size,
                 [&](size_t i) { in_frontier[frontier[i]] = true; });
  }

  void dense2sparse() {
    auto identity = delayed_seq<NodeId>(G.n, [&](size_t i) { return i; });
    pack_into_uninitialized(identity, in_frontier, frontier);
  }

  bool judge(int round) {
    size_t front_out_edges = 0;
    bool both_sparse_dense = false;
    if (!sparse && frontier_size < (1 << 14)) {
      dense2sparse();
      both_sparse_dense = true;
    }
    if (sparse || both_sparse_dense) {
      auto degree_seq = delayed_seq<size_t>(frontier_size, [&](size_t i) {
        NodeId u = frontier[i];
        return G.offsets[u + 1] - G.offsets[u];
      });
      front_out_edges = reduce(degree_seq);
    } else {
      front_out_edges = dense_sample(hash32(round));
    }
    bool sparse_now = ((frontier_size + front_out_edges) < threshold);
    if (!both_sparse_dense && sparse != sparse_now) {
      if (!sparse_now) {
        sparse2dense();
      } else {
        dense2sparse();
      }
    }
    return sparse_now;
  }

 public:
  LDD() = delete;
  LDD(
      const Graph& _G, function<bool(NodeId, NodeId)> _pred =
                           [](NodeId, NodeId) { return true; })
      : G(_G), bag(G.n), pred(_pred) {
    frontier = sequence<NodeId>(G.n);
    in_frontier = sequence<bool>::uninitialized(G.n);
    in_next_frontier = sequence<bool>::uninitialized(G.n);
    sparse = true;
    frontier_size = 0;
    threshold = G.m / 20;
  };

  tuple<sequence<NodeId>, sequence<NodeId>> ldd(double beta,
                                                bool spanning_tree = false) {
    size_t n = G.n;
    sequence<NodeId> label(n, NODE_MAX);
    sequence<NodeId> parent;
    if (spanning_tree) {
      parent = tabulate(n, [&](NodeId i) { return i; });
    }
    sequence<NodeId> perm = sequence<NodeId>(NUM_SAMPLES);
    for (size_t i = 0; i < NUM_SAMPLES; i++) {
      perm[i] = hash32(NUM_SAMPLES + i) % n;
    }
    size_t num_sampled = 0;
    int round = 0;
    frontier_size = 0;
    sparse = true;
    while (frontier_size > 0 || num_sampled < NUM_SAMPLES) {
      round++;
      size_t step_size = floor(exp(round * beta));
      size_t work_size = min(step_size, NUM_SAMPLES - num_sampled);
      size_t num_new_centers = 0;
      if (sparse && work_size > 0) {
        auto centers = filter(perm.cut(num_sampled, num_sampled + work_size),
                              [&](NodeId u) { return label[u] == NODE_MAX; });
        num_new_centers = centers.size();
        parallel_for(0, num_new_centers, [&](size_t i) {
          frontier[frontier_size + i] = centers[i];
          label[centers[i]] = centers[i];
        });
      } else if (work_size > 0) {
        num_new_centers =
            count(perm.cut(num_sampled, num_sampled + work_size), NODE_MAX);
        parallel_for(num_sampled, num_sampled + work_size, [&](NodeId i) {
          NodeId u = perm[i];
          if (label[u] == NODE_MAX) {
            label[u] = u;
            in_frontier[u] = true;
          }
        });
      }
      frontier_size += num_new_centers;
      num_sampled += work_size;
      bool next_sparse = judge(round);
      if (next_sparse) {
        frontier_size = sparse_update(label, parent);
      } else {
        frontier_size = dense_update(label, parent);
      }
      sparse = next_sparse;
    }
    parallel_for(0, G.n, [&](size_t i) {
      if (label[i] == NODE_MAX) {
        label[i] = i;
      }
    });
    return {label, parent};
  }
};
