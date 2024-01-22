#pragma once
#include "connectivity.h"
#include "graph.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "spanning_forest.h"
#include "sparse_table.h"
#include "utils.h"

using namespace std;
using namespace parlay;

template <class Graph, class TagId = typename Graph::NodeId>
class BCC {
 protected:
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr NodeId NODE_MAX = numeric_limits<NodeId>::max();
  static constexpr TagId TAG_MAX = numeric_limits<TagId>::max();
  static constexpr size_t BLOCK_SIZE = 1024;
  static constexpr double beta = 0.2;

  const Graph &G;
  sequence<NodeId> parent;
  sequence<TagId> first;
  sequence<TagId> last;
  sequence<TagId> low;
  sequence<TagId> high;

  sequence<NodeId> euler_tour_tree(const Forest<NodeId> &F) {
    size_t num_trees = F.num_trees;
    size_t m = F.G.m;
    auto edgelist = sequence<pair<NodeId, NodeId>>::uninitialized(m * 2);
    auto perms = sequence<pair<NodeId, TagId>>::uninitialized(m * 2);
    assert(F.G.m == F.G.offsets[F.G.n]);
    parallel_for(0, F.G.n, [&](size_t i) {
      parallel_for(
          F.G.offsets[i], F.G.offsets[i + 1],
          [&](size_t j) {
            edgelist[j * 2] = {F.vertex[i], F.G.edges[j].v};
            edgelist[j * 2 + 1] = {F.G.edges[j].v, F.vertex[i]};
            perms[j * 2] = {F.vertex[i], j * 2};
            perms[j * 2 + 1] = {F.G.edges[j].v, j * 2 + 1};
          },
          BLOCK_SIZE);
    });
    integer_sort_inplace(make_slice(perms),
                         [](const pair<NodeId, TagId> &a) { return a.first; });
    auto first_edge = sequence<TagId>::uninitialized(G.n);
    parallel_for(0, m * 2, [&](size_t i) {
      if (i == 0 || perms[i - 1].first != perms[i].first) {
        first_edge[perms[i].first] = perms[i].second;
      }
    });
    auto link = sequence<TagId>::uninitialized(m * 2);
    parallel_for(0, m * 2, [&](size_t i) {
      if (i + 1 < m * 2 && perms[i].first == perms[i + 1].first) {
        link[perms[i].second ^ 1] = perms[i + 1].second;
      } else {
        link[perms[i].second ^ 1] = first_edge[perms[i].first];
      }
    });

    auto samples_offsets = sequence<TagId>::uninitialized(num_trees + 1);
    parallel_for(0, num_trees, [&](size_t i) {
      size_t tree_size =
          (F.G.offsets[F.offsets[i + 1]] - F.G.offsets[F.offsets[i]]) + 1;
      size_t edges_size = 2 * (tree_size - 1);
      samples_offsets[i] = sqrt(edges_size);
    });
    samples_offsets[num_trees] = 0;
    size_t num_samples = scan_inplace(make_slice(samples_offsets));
    auto samples = sequence<TagId>::uninitialized(num_samples);
    auto skip_to = sequence<pair<TagId, TagId>>::uninitialized(num_samples);
    sequence<TagId> idx(m * 2, TAG_MAX);
    parallel_for(0, num_trees, [&](size_t i) {
      size_t tree_size =
          (F.G.offsets[F.offsets[i + 1]] - F.G.offsets[F.offsets[i]]) + 1;
      size_t edges_size = 2 * (tree_size - 1);
      for (size_t j = samples_offsets[i]; j < samples_offsets[i + 1]; j++) {
        if (j == samples_offsets[i]) {
          samples[j] = 2 * (TagId)F.G.offsets[F.offsets[i]];
        } else {
          uint64_t seed = hash64(j);
          TagId pos = seed % edges_size;
          while (idx[pos + 2 * (TagId)F.G.offsets[F.offsets[i]]] != TAG_MAX) {
            pos = (pos + 1) % edges_size;
          }
          samples[j] = pos + 2 * (TagId)F.G.offsets[F.offsets[i]];
        }
        idx[samples[j]] = j;
      }
    });

    parallel_for(0, num_trees, [&](size_t i) {
      parallel_for(samples_offsets[i], samples_offsets[i + 1], [&](size_t j) {
        TagId node = samples[j];
        skip_to[j].second = 0;
        do {
          node = link[node];
          skip_to[j].second++;
        } while (idx[node] == TAG_MAX);
        skip_to[j].first = idx[node];
      });
    });

    auto offsets = sequence<TagId>::uninitialized(num_samples);
    parallel_for(0, num_trees, [&](size_t i) {
      size_t sum = 0;
      TagId cur_idx = samples_offsets[i];
      for (size_t j = samples_offsets[i]; j < samples_offsets[i + 1]; j++) {
        offsets[cur_idx] = sum;
        sum += skip_to[cur_idx].second;
        cur_idx = skip_to[cur_idx].first;
      }
    });

    auto sizes = sequence<TagId>::uninitialized(num_trees + 1);
    parallel_for(0, num_trees, [&](size_t i) {
      sizes[i] = 2 * (TagId)(F.G.offsets[F.offsets[i + 1]] -
                             F.G.offsets[F.offsets[i]]) +
                 1;
    });
    sizes[num_trees] = 0;
    scan_inplace(make_slice(sizes));
    auto order = sequence<NodeId>::uninitialized(sizes[num_trees]);

    parallel_for(0, num_trees, [&](size_t i) {
      parallel_for(samples_offsets[i], samples_offsets[i + 1], [&](size_t j) {
        TagId node = samples[j];
        TagId id = idx[node];
        TagId cur_offsets = offsets[id];
        do {
          order[sizes[i] + cur_offsets] = edgelist[node].first;
          cur_offsets++;
          node = link[node];
        } while (idx[node] == TAG_MAX);
      });
      if (samples_offsets[i] != samples_offsets[i + 1]) {
        order[sizes[i + 1] - 1] = edgelist[samples[samples_offsets[i]]].first;
      } else {
        order[sizes[i + 1] - 1] = F.vertex[F.offsets[i]];
      }
    });
    return order;
  }

  void tagging(const Forest<NodeId> &F, const sequence<NodeId> &order) {
    size_t n = G.n;
    parallel_for(0, n, [&](size_t i) {
      first[i] = TAG_MAX;
      last[i] = 0;
      parent[i] = i;
    });
    parallel_for(0, order.size(), [&](TagId i) {
      NodeId v = order[i];
      write_min(&first[v], i);
      write_max(&last[v], i);
    });
    parallel_for(0, F.G.n, [&](size_t i) {
      NodeId u = F.vertex[i];
      parallel_for(
          F.G.offsets[i], F.G.offsets[i + 1],
          [&](size_t j) {
            NodeId v = F.G.edges[j].v;
            if (first[u] < first[v]) {
              parent[v] = u;
            } else {
              parent[u] = v;
            }
          },
          BLOCK_SIZE);
    });

    auto w = tabulate(first.size(),
                      [&](size_t i) { return make_pair(first[i], first[i]); });
    parallel_for(0, n, [&](size_t i) {
      parallel_for(
          G.offsets[i], G.offsets[i + 1],
          [&](size_t j) {
            NodeId u = i, v = G.edges[j].v;
            if (u < v && parent[u] != v && parent[v] != u) {
              if (first[u] < first[v]) {
                write_min(&w[v].first, first[u]);
                write_max(&w[u].second, first[v]);
              } else {
                write_min(&w[u].first, first[v]);
                write_max(&w[v].second, first[u]);
              }
            }
          },
          BLOCK_SIZE);
    });

    auto seq = tabulate(order.size(), [&](size_t i) { return w[order[i]]; });
    sparse_table st(seq, minmaxm<TagId>());

    parallel_for(0, n, [&](size_t i) {
      std::tie(low[i], high[i]) = st.query(first[i], last[i] + 1);
    });
  }

  void get_component_head(const sequence<NodeId> &label) {
    size_t n = G.n;
    auto component_head = sequence<NodeId>(n, NODE_MAX);
    parallel_for(0, n, [&](NodeId i) {
      NodeId p = parent[i];
      if (label[p] != label[i]) {
        component_head[label[i]] = p;
      }
    });
  }

 public:
  BCC() = delete;
  BCC(const Graph &_G) : G(_G) {
    parent = sequence<NodeId>::uninitialized(G.n);
    first = sequence<TagId>::uninitialized(G.n);
    last = sequence<TagId>::uninitialized(G.n);
    low = sequence<TagId>::uninitialized(G.n);
    high = sequence<TagId>::uninitialized(G.n);
  }

  sequence<NodeId> biconnectivity() {
    Forest F = spanning_forest(G, beta);
    sequence<NodeId> order = euler_tour_tree(F);
    tagging(F, order);
    auto critical = [&](NodeId u, NodeId v) {
      if (first[u] <= low[v] && last[u] >= high[v]) {
        return true;
      }
      return false;
    };
    auto backward = [&](NodeId u, NodeId v) {
      return first[u] <= first[v] && last[u] >= first[v];
    };
    function<bool(NodeId, NodeId)> pred = [&](NodeId u, NodeId v) -> bool {
      if (parent[v] == u) {
        if (!critical(u, v)) {
          return true;
        }
      } else if (parent[u] == v) {
        if (!critical(v, u)) {
          return true;
        }
      } else {
        if (!backward(u, v) && !backward(v, u)) {
          return true;
        }
      }
      return false;
    };

    auto label = get<0>(connectivity(G, beta, pred));
    get_component_head(label);
    return label;
  }

  void get_num_bcc(const sequence<NodeId> &label) {
    auto cc_label = get<0>(connectivity(G, beta));
    auto unique_cc_label = remove_duplicates_ordered(cc_label, less<NodeId>());
    auto unique_bcc_label = remove_duplicates_ordered(label, less<NodeId>());
    auto sorted_label = sort(label);
    size_t n = sorted_label.size();
    auto unique_index = pack_index(delayed_seq<bool>(n + 1, [&](size_t i) {
      return i == 0 || i == n || sorted_label[i] != sorted_label[i - 1];
    }));
    auto largest_bcc = reduce(
        delayed_seq<size_t>(
            unique_index.size() - 1,
            [&](size_t i) { return unique_index[i + 1] - unique_index[i]; }),
        maxm<size_t>());
    printf("#CC: %zu\n", unique_cc_label.size());
    printf("#BCC: %zu\n", unique_bcc_label.size() - unique_cc_label.size());
    printf("Largest_BCC: %zu\n", largest_bcc);
    ofstream ofs("fast-bcc.tsv", ios_base::app);
    ofs << unique_cc_label.size() << '\t'
        << unique_bcc_label.size() - unique_cc_label.size() << '\t';
    ofs.close();
  }
};