#include <set>

#include "graph.h"
#include "hashbag.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

using namespace std;
using namespace parlay;

template <typename E, typename EV>
inline pair<E, bool> fetch_and_add_bounded(E *a, EV b, E k) {
  volatile E newV, oldV;
  bool c = false;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (oldV > k && !(c = atomic_compare_and_swap(a, oldV, newV)));
  return make_pair(oldV, c);
}

template <class Graph>
class KCore {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  // tunable parameters
  static constexpr bool enable_sampling = true;
  static constexpr bool enable_local_queue = true;
  static constexpr uint32_t log2_single_buckets = 3;
  static constexpr uint32_t num_intermediate_buckets = 6;
  static constexpr uint32_t sample_threshold = 2000;
  static constexpr size_t BLOCK_SIZE = 1024;
  static constexpr double init_reduce_ratio = 0.1;
  static constexpr uint32_t log2_error_factor = 32;
  static constexpr uint32_t local_queue_size = 128;
  static constexpr double bias_factor = 0.5;
  static constexpr double error_rate_tolerance = 0.0000000001;

  // other parameters
  static constexpr uint32_t num_single_buckets = 1 << log2_single_buckets;
  static constexpr uint32_t bucket_mask = num_single_buckets - 1;
  static constexpr uint32_t stride = num_single_buckets
                                     << num_intermediate_buckets;
  static constexpr uint32_t exp_hits =
      log2_error_factor / (init_reduce_ratio * init_reduce_ratio);

  const Graph &G;
  sequence<hashbag<NodeId>> buckets;
  sequence<NodeId> frontier;
  sequence<NodeId> coreness;
  sequence<bool> alive;
  sequence<bool> sample_mode;
  sequence<Sampler> samplers;
  hashbag<NodeId> counting_bag;

 public:
  KCore() = delete;
  KCore(const Graph &_G) : G(_G), counting_bag(G.n) {
    size_t n = G.n;
    buckets = sequence<hashbag<NodeId>>(
        num_single_buckets + num_intermediate_buckets, hashbag<NodeId>(n));
    frontier = sequence<NodeId>::uninitialized(n);
    coreness = sequence<NodeId>::uninitialized(n);
    alive = sequence<bool>::uninitialized(n);
    sample_mode = sequence<bool>::uninitialized(n);
    samplers = sequence<Sampler>::uninitialized(n);
  }

  // make sure coreness[v] is the accurate remaining degree before calling the
  // function
  void set_sampler(NodeId v, NodeId k) {
    if (coreness[v] * init_reduce_ratio >= sample_threshold &&
        k < coreness[v] * init_reduce_ratio * bias_factor &&
        exp_hits < (coreness[v] - k) * (1 - init_reduce_ratio)) {
      sample_mode[v] = true;
      // size_t modified_exp_hits = std::min((size_t)((coreness[v] - k) / 16),
      // (size_t)exp_hits);
      size_t modified_exp_hits = log2_error_factor *
                                 parlay::log2_up((int)((coreness[v] - k))) *
                                 parlay::log2_up((int)((coreness[v] - k))) / 2;
      // cout << "modified_exp_hits: " << modified_exp_hits << endl;
      double sample_rate =
          modified_exp_hits / ((1 - init_reduce_ratio) * coreness[v]);
      samplers[v].reset(modified_exp_hits, sample_rate);
    } else {
      sample_mode[v] = false;
    }
  }

  // add a vertex u with degree d to the corresponding bucket
  void add_to_bucket(NodeId u, NodeId d, NodeId base_k) {
    // degree out of range
    if (d < base_k || d > (base_k | (stride - 1))) {
      return;
    }
    if (d < base_k + num_single_buckets) {
      // add to single bucket
      buckets[d & bucket_mask].insert(u);
    } else {
      // add to intermediate buckets
      NodeId diff_bit = 63 - __builtin_clzll(d ^ base_k);
      buckets[diff_bit - log2_single_buckets + num_single_buckets].insert(u);
    }
  }

  inline void move_bucket(NodeId u, NodeId d, NodeId base_k) {
    // degree out of range
    if (d < base_k || d > (base_k | (stride - 1))) {
      return;
    }
    // add to single bucket
    if (d < base_k + num_single_buckets) {
      buckets[d & bucket_mask].insert(u);
      return;
    }
    // add to intermediate buckets
    NodeId diff_bit = 63 - __builtin_clzll(d ^ base_k);
    NodeId previous_diff_bit = 63 - __builtin_clzll((d + 1) ^ base_k);
    if (diff_bit != previous_diff_bit) {
      buckets[diff_bit - log2_single_buckets + num_single_buckets].insert(u);
    }
  }

  void sample_vertex(NodeId u, NodeId v, bool &counting_flag) {
    uint32_t hash_v = hash32(u * G.n + v);
    bool callback = false;
    samplers[v].sample(hash_v, callback);
    if (callback) {
      // the degree has been reduced by reduce_ratio
      if (counting_flag == false) {
        counting_flag = true;
      }
      counting_bag.insert(v);
    }
  }

  void fetch_and_add_vertex(NodeId v, NodeId base_k, NodeId k) {
    auto [id, succeed] = fetch_and_add_bounded(&coreness[v], -1, k);
    id--;
    if (succeed) {
      move_bucket(v, id, base_k);
    }
  }

  void count_alive_neighbors(NodeId u) {
    coreness[u] = parlay::count_if(G.edges.cut(G.offsets[u], G.offsets[u + 1]),
                                   [&](auto &es) { return alive[es.v]; });
  }

  void count_vertex(NodeId u, NodeId k, NodeId base_k) {
    NodeId was = coreness[u];
    count_alive_neighbors(u);
    if (coreness[u] < k) {
      NodeId alive_last_round = parlay::count_if(
          G.edges.cut(G.offsets[u], G.offsets[u + 1]),
          [&](auto &es) { return alive[es.v] || coreness[es.v] == k; });
      if (alive_last_round >= k) {
        // deg[u] is reduced to k in this round
        coreness[u] = k;
        buckets[k & bucket_mask].insert(u);
      } else {
        // deg[u] is reduced to k in previous rounds. Fatal error.
        printf("coreness[%u]: %u, was: %u, k: %u, alive_last_round: %u\n", u,
               coreness[u], was, k, alive_last_round);
        assert(coreness[u] >= k);
      }
    } else {
      add_to_bucket(u, coreness[u], base_k);
    }
    // reset the sampler after counting the actual degree
    set_sampler(u, k);
  }

  void count_vertex_wo_bucketing(NodeId u, NodeId k) {
    NodeId was = coreness[u];
    count_alive_neighbors(u);
    if (coreness[u] < k) {
      NodeId alive_last_round = parlay::count_if(
          G.edges.cut(G.offsets[u], G.offsets[u + 1]),
          [&](auto &es) { return alive[es.v] || coreness[es.v] == k; });
      if (alive_last_round >= k) {
        // deg[u] is reduced to k in this round
        coreness[u] = k;
        buckets[k & bucket_mask].insert(u);
      } else {
        // deg[u] is reduced to k in previous rounds. Fatal error.
        printf("coreness[%u]: %u, was: %u, k: %u, alive_last_round: %u\n", u,
               coreness[u], was, k, alive_last_round);
        assert(coreness[u] >= k);
      }
    }
    // reset the sampler after counting the actual degree
    set_sampler(u, k);
  }

  void map_neighbors_parallel(NodeId u, NodeId base_k, NodeId k,
                              bool &counting_flag) {
    parallel_for(G.offsets[u], G.offsets[u + 1], [&](size_t es) {
      auto v = G.edges[es].v;
      if (coreness[v] > k) {
        if (enable_sampling && sample_mode[v]) {
          sample_vertex(u, v, counting_flag);
        } else {
          fetch_and_add_vertex(v, base_k, k);
        }
      }
    });
  }

  void map_neighbors_sequential(NodeId u, NodeId base_k, NodeId k,
                                bool &counting_flag, NodeId *local_queue,
                                size_t &rear) {
    for (EdgeId i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      NodeId v = G.edges[i].v;
      if (coreness[v] > k) {
        if (enable_sampling && sample_mode[v]) {
          sample_vertex(u, v, counting_flag);
        } else {
          auto [id, succeed] = fetch_and_add_bounded(&coreness[v], -1, k);
          id--;
          if (succeed) {
            if (enable_local_queue && id == k && rear < local_queue_size) {
              local_queue[rear++] = v;
            } else {
              move_bucket(v, id, base_k);
            }
          }
        }
      }
    }
  }

  void map_neighbors_sequentially_wo_bucketing(NodeId u, NodeId k,
                                               bool &counting_flag,
                                               NodeId *local_queue,
                                               size_t &rear) {
    for (EdgeId i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
      NodeId v = G.edges[i].v;
      if (coreness[v] > k) {
        if (enable_sampling && sample_mode[v]) {
          sample_vertex(u, v, counting_flag);
        } else {
          auto [id, succeed] = fetch_and_add_bounded(&coreness[v], -1, k);
          id--;
          if (succeed && id == k) {
            if (enable_local_queue && id == k && rear < local_queue_size) {
              local_queue[rear++] = v;
            } else {
              // insert to the first
              buckets[0].insert(v);
            }
          }
        }
      }
    }
  }

  void map_neighbors_parallel_wo_bucketing(NodeId u, NodeId k,
                                           bool &counting_flag) {
    parallel_for(G.offsets[u], G.offsets[u + 1], [&](size_t es) {
      auto v = G.edges[es].v;
      if (coreness[v] > k) {
        if (enable_sampling && sample_mode[v]) {
          sample_vertex(u, v, counting_flag);
        } else {
          auto [id, succeed] = fetch_and_add_bounded(&coreness[v], -1, k);
          id--;
          if (succeed && id == k) {
            // insert to the first
            buckets[0].insert(v);
          }
        }
      }
    });
  }

  double check_sample_security(NodeId v, size_t k, double sample_rate) {
    double error_probability_bound = 0;
    if (coreness[v] * init_reduce_ratio * bias_factor < k) {
      return 1;
    } else {
      int n_star = coreness[v] - k;
      size_t num_hits = std::max((uint32_t)1, samplers[v].get_num_hits());
      error_probability_bound =
          std::exp(-1.0 * n_star * sample_rate + 2 * num_hits -
                   1.0 * num_hits * num_hits / n_star / sample_rate);
      return error_probability_bound;
    }
  }

  sequence<NodeId> kcore() {
    size_t n = G.n;
    size_t bucketing_pt = 16;
    auto remaining_vertices = parlay::sequence<NodeId>::uninitialized(n);
    size_t avg_deg = G.m / n;
    parallel_for(0, n, [&](size_t i) { remaining_vertices[i] = i; });
    bool contains_sampling_nodes = false;
    // init
    parallel_for(0, n, [&](size_t i) {
      coreness[i] = G.offsets[i + 1] - G.offsets[i];
      if (enable_sampling &&
          coreness[i] * init_reduce_ratio >= sample_threshold) {
        contains_sampling_nodes = true;
      }
      alive[i] = true;
    });
    NodeId max_core = 0;

    if (contains_sampling_nodes) {
      parallel_for(0, n, [&](size_t i) { set_sampler(i, 0); });
    }

    internal::timer t_insert("insert", false);
    internal::timer t_dump("dump", false);
    internal::timer t_push("push", false);

    internal::timer t_pack("pack", false);
    internal::timer t_add("add", false);
    internal::timer t_check_n_count("check_n_count", false);
    internal::timer t_rho_round("round", false);
    size_t num_rho = 0;

    // process from 0 to 16 using a single bucket
    // for (NodeId k = 0; k < bucketing_pt; k++) {
    size_t k = 0;
    if (avg_deg < bucketing_pt) {
      while (k < bucketing_pt) {
        // size_t sub_rho = 0;
        if (remaining_vertices.size() == 0) {
          break;
        }
        t_insert.start();
        parallel_for(0, remaining_vertices.size(), [&](size_t i) {
          // add to the first bucket
          if (coreness[remaining_vertices[i]] == k) {
            buckets[0].insert(remaining_vertices[i]);
          }
          if (contains_sampling_nodes && sample_mode[remaining_vertices[i]]) {
            // prob check
            t_check_n_count.start();
            double error_rate = check_sample_security(
                remaining_vertices[i], k + stride,
                samplers[remaining_vertices[i]].get_exp_hits() /
                    ((1 - init_reduce_ratio) *
                     coreness[remaining_vertices[i]]));
            if (error_rate >= error_rate_tolerance) {
              t_check_n_count.start();
              count_vertex(remaining_vertices[i], k, k);
            }
            t_check_n_count.stop();
          }
        });
        t_insert.stop();
        // pack bucket 0
        t_pack.start();
        auto size = buckets[0].pack_into(make_slice(frontier));
        t_pack.stop();
        t_push.start();
        while (size) {
          // sub_rho++;
          num_rho++;
          bool counting_flag = false;
          parallel_for(
              0, size,
              [&](size_t j) {
                auto f = frontier[j];
                alive[f] = false;
                if (max_core < coreness[f]) {
                  max_core = coreness[f];
                }
                if (enable_local_queue) {
                  NodeId local_queue[local_queue_size];
                  size_t front = 0, rear = 0;
                  local_queue[rear++] = f;
                  while (front < rear) {
                    NodeId u = local_queue[front++];
                    alive[u] = false;
                    size_t deg = G.offsets[u + 1] - G.offsets[u];
                    if (deg < BLOCK_SIZE) {
                      // sequentially insert
                      map_neighbors_sequentially_wo_bucketing(
                          u, k, counting_flag, local_queue, rear);
                    } else {
                      // parallel insert
                      map_neighbors_parallel_wo_bucketing(u, k, counting_flag);
                    }
                  }
                } else {
                  map_neighbors_parallel_wo_bucketing(f, k, counting_flag);
                }
              },
              1);
          if (counting_flag) {
            t_check_n_count.start();
            auto counting_vertices = counting_bag.pack();
            parallel_for(0, counting_vertices.size(), [&](size_t j) {
              NodeId u = counting_vertices[j];
              count_vertex_wo_bucketing(u, k);
            });
            t_check_n_count.stop();
          }
          size = buckets[0].pack_into(make_slice(frontier));
        }
        t_push.stop();
        remaining_vertices = parlay::filter(remaining_vertices,
                                            [&](NodeId v) { return alive[v]; });
        k++;
      }
    }

    // remaining vertices using hierrachical buckets
    if (remaining_vertices.size() > 0) {
      for (NodeId base_k = 0;; base_k += stride) {
        if (remaining_vertices.size() == 0) {
          break;
        }
        t_insert.start();
        parallel_for(0, remaining_vertices.size(), [&](size_t i) {
          add_to_bucket(remaining_vertices[i], coreness[remaining_vertices[i]],
                        base_k);
          if (contains_sampling_nodes && sample_mode[remaining_vertices[i]]) {
            // prob check
            t_check_n_count.start();
            double error_rate = check_sample_security(
                remaining_vertices[i], base_k + stride,
                samplers[remaining_vertices[i]].get_exp_hits() /
                    ((1 - init_reduce_ratio) *
                     coreness[remaining_vertices[i]]));
            if (error_rate >= error_rate_tolerance) {
              count_vertex(remaining_vertices[i], base_k, base_k);
            }
            t_check_n_count.stop();
          }
        });
        t_insert.stop();
        NodeId offset_k = 0;
        for (NodeId k = base_k; k < base_k + stride; k++) {
          // size_t sub_rho = 0;

          t_dump.start();
          if (base_k != k) {
            for (int i = num_intermediate_buckets - 1; i >= 0; i--) {
              uint32_t mask = (num_single_buckets << i) - 1;
              if ((k & mask) == 0) {
                offset_k += num_single_buckets;
                t_pack.start();
                // printf("packing\n");
                auto size = buckets[num_single_buckets + i].pack_into(
                    make_slice(frontier));
                // printf("dumping %zu vertices\n", size);
                t_pack.stop();
                t_add.start();
                parallel_for(0, size, [&](size_t j) {
                  NodeId u = frontier[j];
                  add_to_bucket(u, coreness[u], base_k + offset_k);
                });
                t_add.stop();
                // for(size_t j = 0; j < num_intermediate_buckets; j++) {
                // printf("buckets %zu: \n", j);
                // buckets[j + num_single_buckets].print();
                //}
                break;
              }
            }
          }
          t_dump.stop();
          // printf("k: %u, base_k: %u\n", k, base_k + offset_k);
          t_push.start();
          auto size = buckets[k & bucket_mask].pack_into_pred(
              make_slice(frontier), [&](NodeId v) { return coreness[v] == k; });
          while (size) {
            // sub_rho++;
            num_rho++;
            // t_rho_round.start();
            bool counting_flag = false;
            parallel_for(
                0, size,
                [&](size_t j) {
                  auto f = frontier[j];
                  alive[f] = false;
                  if (max_core < coreness[f]) {
                    max_core = coreness[f];
                  }
                  if (enable_local_queue) {
                    NodeId local_queue[local_queue_size];
                    size_t front = 0, rear = 0;
                    local_queue[rear++] = f;
                    while (front < rear) {
                      NodeId u = local_queue[front++];
                      alive[u] = false;
                      size_t deg = G.offsets[u + 1] - G.offsets[u];
                      if (deg < BLOCK_SIZE) {
                        map_neighbors_sequential(u, base_k + offset_k, k,
                                                 counting_flag, local_queue,
                                                 rear);
                      } else {
                        map_neighbors_parallel(u, base_k + offset_k, k,
                                               counting_flag);
                      }
                    }
                  } else {
                    map_neighbors_parallel(f, base_k + offset_k, k,
                                           counting_flag);
                  }
                },
                1);
            if (counting_flag) {
              auto counting_vertices = counting_bag.pack();
              t_check_n_count.start();
              parallel_for(0, counting_vertices.size(), [&](size_t j) {
                NodeId u = counting_vertices[j];
                t_check_n_count.start();
                count_vertex(u, k, base_k + offset_k);
                t_check_n_count.stop();
              });
              t_check_n_count.stop();
            }
            size = buckets[k & bucket_mask].pack_into_pred(
                make_slice(frontier),
                [&](NodeId v) { return coreness[v] == k; });
            // t_rho_round.next("round");
          }
          t_push.stop();
          // if (k >= 8) {
          //   cout << sub_rho << endl;
          // }
        }
        t_pack.start();
        remaining_vertices = parlay::filter(remaining_vertices,
                                            [&](NodeId v) { return alive[v]; });
        t_pack.stop();
      }
    }
    printf("coreness: %u\n", max_core);
    t_insert.total();
    t_dump.total();
    t_push.total();
    t_pack.total();
    t_add.total();
    t_check_n_count.total();
    cout << "rho: " << num_rho << endl;
    return coreness;
  }
};