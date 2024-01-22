#pragma once
#include <climits>
#include <queue>

#include "graph.h"
#include "hashbag.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "utils.h"


using namespace parlay;
using namespace std;


template <class Graph>
class REACH {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr size_t LOCAL_QUEUE_SIZE = 512;
  static constexpr size_t SPARSE_TH = 128;
  static constexpr size_t BLOCK_SIZE = 1024;


  hashbag<NodeId> bag;
  sequence<NodeId> frontier;
  sequence<bool> dense_frontier;
  NodeId n_frontier;

 public:
  size_t num_round;
  REACH() = delete;
  REACH(size_t n) : bag(n) {
    frontier=sequence<NodeId>::uninitialized(n);
    dense_frontier = sequence<bool>::uninitialized(n);
  }


  size_t sparse_visit(const Graph& G, sequence<bool>& vist) {
    parallel_for(0, n_frontier,[&](size_t i) {
        NodeId f = frontier[i];
        size_t f_visit = 0;
        size_t deg_f = G.offsets[f + 1] - G.offsets[f];
        if ((deg_f < LOCAL_QUEUE_SIZE) && (deg_f > 0) &&
            (f_visit < LOCAL_QUEUE_SIZE)) {
          NodeId Q[LOCAL_QUEUE_SIZE];
          size_t head = 0, tail = 0;
          Q[tail++] = f;
          while (head < tail && tail < LOCAL_QUEUE_SIZE) {
            NodeId u = Q[head++];
            size_t deg_u = G.offsets[u + 1] - G.offsets[u];
            if (deg_u > LOCAL_QUEUE_SIZE) {
              bag.insert(u);
              break;
            }
            for (size_t j = G.offsets[u]; j < G.offsets[u + 1]; j++) {
              NodeId v = G.edges[j].v;
              if (!vist[v] && compare_and_swap(&vist[v], false, true)) {
                f_visit++;
                if (f_visit < LOCAL_QUEUE_SIZE) {
                  Q[tail++] = v;
                } else {
                  bag.insert(v);
                }
              }
            }
          }
          for (size_t j = head; j < tail; j++) {
            bag.insert(Q[j]);
          }
        } else if (deg_f > 0) {
          parallel_for(G.offsets[f], G.offsets[f + 1],[&](size_t j) {
            if (j>=G.edges.size()){
              printf("offsets[f]:%lu offsets[f+1]&:%lu j: %zu, f:%u\n",G.offsets[f], G.offsets[f+1], j, f);
            }
            assert(j<G.edges.size());
            NodeId v = G.edges[j].v;
            if (!vist[v] && compare_and_swap(&vist[v], false, true)){bag.insert(v);}
          },BLOCK_SIZE);
        }
      },1);
    return bag.pack_into(make_slice(frontier));
  }

  size_t dense_visit(const Graph& GT, sequence<bool>& vist) {
    parallel_for(0, GT.n, [&](size_t i) {
      if (vist[i] == false) {
        dense_frontier[i]=false;
        for (size_t j = GT.offsets[i]; j < GT.offsets[i + 1]; j++) {
          NodeId ngb_node = GT.edges[j].v;
          if (vist[ngb_node]) {
            vist[i] = true;
            dense_frontier[i] = true;
            break;}}
      }else{
        dense_frontier[i]=false;
      }
    },BLOCK_SIZE);
    return parlay::count(make_slice(dense_frontier), true);
  }

 void reach(NodeId source, const Graph& G, const Graph& GT, parlay::sequence<bool>&visit) {
    parlay::parallel_for(0, G.n, [&](size_t i){visit[i]=false;});
    n_frontier = 0;
    frontier[n_frontier++] = source;
    visit[source] = true;
    bool is_sparse = true;       // whether is sparse edge map
    num_round = 0;
    while (n_frontier > 0) {
      num_round++;
      if (is_sparse) {
        auto n_edges = parlay::reduce(parlay::delayed_map(frontier.cut(0, n_frontier), 
                        [&] (NodeId i) {return G.offsets[i+1]-G.offsets[i];}));
        if ((n_frontier + n_edges) > G.m/10) {
          parallel_for(0, G.n, [&](NodeId i){dense_frontier[i]=false;});
          parlay::for_each(frontier.cut(0, n_frontier),[&](NodeId v){dense_frontier[v]=true;});
          is_sparse=false;
        } else is_sparse =true;
      } else {
        if (n_frontier > G.n/20) is_sparse = false;
        else {
          // frontier = parlay::pack_index<NodeId>(dense_frontier);
          auto identity = [](size_t i) -> NodeId { return static_cast<NodeId>(i); };
          parlay::pack_into_uninitialized(parlay::delayed_tabulate(G.n,identity),dense_frontier, frontier);
          is_sparse = true;}
      }
      n_frontier = (is_sparse)? sparse_visit(G, visit): dense_visit(GT, visit);
    }
  }
  void verifier(const Graph& G, NodeId s, const sequence<bool> &act_dist) {
    internal::timer t;
    size_t n = G.n;
    sequence<bool> exp_dist(n, false);
    exp_dist[s] = true;
    queue<NodeId> q;
    q.push(s);
    while (!q.empty()) {
      NodeId u = q.front();
      q.pop();
      for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
        NodeId v = G.edges[i].v;
        if (!exp_dist[v]) {
          exp_dist[v] =true;
          q.push(v);
        }
      }
    }
    t.stop();
    printf("Sequential BFS running time: %f\n", t.total_time());

    assert(exp_dist.size() == act_dist.size());
    for (size_t i = 0; i < n; i++) {
      if (exp_dist[i] != act_dist[i]) {
        printf("exp_dist[%zu]: %u whlie act_dist[%zu]: %u\n", i, exp_dist[i], i,
              act_dist[i]);
      }
      assert(exp_dist[i] == act_dist[i]);
    }
    printf("Passed!\n");
  }
};
