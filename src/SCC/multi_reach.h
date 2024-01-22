#pragma once
#include <climits>
#include <queue>

#include "graph.h"
#include "hashbag.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "utils.h"

#include "resizable_table.h"

#include "reach.h"



using namespace parlay;
using namespace std;


template <class Graph, typename LabelTy>
class MULTI_REACH {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;

  static constexpr size_t LOCAL_QUEUE_SIZE = 512;

  hashbag<NodeId> bag;
  sequence<NodeId> frontier;
  sequence<bool> bits;
  NodeId n_frontier;
  size_t n,m;  // n is the remaining vertices, m is the total vertices


 public:
  size_t r; // number of rounds
  MULTI_REACH() = delete;
  MULTI_REACH(size_t _n, size_t _m) : bag(_n), n(_n), m(_m) {
    frontier=sequence<NodeId>::uninitialized(n);
    bits=sequence<bool>(m,false);
  }
  bool multi_reach(const Graph& G, sequence<LabelTy>& label,sequence<NodeId>& sources,
                    hash_table<NodeId, NodeId>& table, LabelTy l){
    n_frontier = sources.size();
    assert(table.m>n_frontier);
    parallel_for(0, n_frontier, [&](size_t i){
      frontier[i]=sources[i];
      table.insert(frontier[i], l+i);
    });
    auto propagate = [&](NodeId u, NodeId v){
      bool label_changed = false;
      size_t num_pairs = 0;
      iter_k<NodeId> u_iter={u,0,0};
      if (table.init_iter(u_iter)){
        while(true){
          NodeId u_label=table.get(u_iter);
          label_changed |= table.insert(v, u_label);
          if (table.overfull){return std::pair(false,(size_t)LOCAL_QUEUE_SIZE);}
          num_pairs++;
          if (!table.has_next(u_iter)){break;}
        }
      }
      return std::pair(label_changed, num_pairs);
    };
    r=0;
    parlay::internal::timer t("multi_search", false);
    while(n_frontier>0){
      r++;
      parallel_for(0, n_frontier, [&](NodeId i){
        NodeId f = frontier[i];
        EdgeId deg_f = G.offsets[f+1]-G.offsets[f];
        if (deg_f>LOCAL_QUEUE_SIZE){
          parlay::for_each(G.edges.cut(G.offsets[f], G.offsets[f+1]),[&](auto e){
            if (label[e.v] == label[f]){
              bool suc=(r==1)?table.insert(e.v,l+i):propagate(f,e.v).first;
              if (suc&&compare_and_swap(&bits[e.v],false,true)){bag.insert(e.v);}
            }
          });
        }else{
          NodeId Q[LOCAL_QUEUE_SIZE];
          NodeId head=0;
          NodeId tail=0;
          Q[tail++]=f;
          size_t f_visit=0;
          while (head < tail && f_visit<LOCAL_QUEUE_SIZE){
            NodeId u = Q[head];
            EdgeId deg_u = G.offsets[u+1]-G.offsets[u];
            if (deg_u>LOCAL_QUEUE_SIZE){break;}
            head++;
            for (size_t j = G.offsets[u]; j<G.offsets[u+1]; j++){
              auto e = G.edges[j];
              if (label[u]==label[e.v]){
                bool suc=false;
                if(r==1){
                  suc=table.insert(e.v, l+i);
                  f_visit++;
                }else{
                  auto ans = propagate(u,e.v);
                  suc=ans.first;
                  f_visit+=ans.second;
                }
                if (suc){
                  if (f_visit<LOCAL_QUEUE_SIZE){
                    Q[tail++]=e.v;
                  }else{
                    if (compare_and_swap(&bits[e.v],false,true)){bag.insert(e.v);}
                  }               
                }
              }
            }
          }
          for(NodeId i = head; i<tail; i++){
            if (compare_and_swap(&bits[Q[i]],false,true))
              bag.insert(Q[i]);
          }
        }
      },1);
      if (table.overfull){return false;}
      n_frontier = bag.pack_into(frontier);
      parlay::for_each(frontier.cut(0, n_frontier),[&](NodeId v){bits[v]=false;});
      t.next("subround time");
    }
    return true;
  }

  int multi_reach_safe(const Graph& G, sequence<LabelTy>& label,sequence<NodeId>& sources,
                    hash_table<NodeId, NodeId>& table, LabelTy l){
    int round=0;
    
    while (true){
      bool succeed = multi_reach(G,label, sources,table,l);
      if (succeed){break;}
      cout << "trigger table resize" << endl;
      parallel_for(0, G.n, [&](size_t i) { bits[i] = false; });
      table.double_size();
      round++;
    }
    return round;
  }

  bool verifier(Graph& G,sequence<LabelTy>& labels, sequence<NodeId>& sources,
      hash_table<NodeId, NodeId>& table,  LabelTy l){
    size_t n = G.n;
    auto pairs = table.pack();
    parlay::sort_inplace(pairs, [&](auto a, auto b){
      return (a.second != b.second)? a.second<b.second:
        a.first<b.first;
    });
    size_t label_i =0;
    parlay::sequence<bool> exp_dist(n, false);
    parlay::sequence<bool> ans_dist(n, false);
    for (size_t si = 0; si < sources.size(); si++){
      NodeId s = sources[si];
      parlay::parallel_for(0, n,[&](size_t i){exp_dist[i]=false;});
      exp_dist[s] = true;
      queue<NodeId> q;
      q.push(s);
      while (!q.empty()) {
        NodeId u = q.front();
        q.pop();
        for (size_t i = G.offsets[u]; i < G.offsets[u + 1]; i++) {
          NodeId v = G.edges[i].v;
          if (labels[u]!=labels[v])continue;
          if (!exp_dist[v]) {
            exp_dist[v] =true;
            q.push(v);
          }
        }
      }
      parlay::parallel_for(0, n, [&](size_t i){ans_dist[i]=false;});
      for (; label_i<pairs.size(); label_i++){
        if (pairs[label_i].second!= si+l){
          break;
        }
        ans_dist[pairs[label_i].first]=true;
      }
      bool pass = true;
      parlay::parallel_for(0, n,[&](size_t i){
        if (ans_dist[i] != exp_dist[i]){
          pass &=false;
        }
      });
      if (!pass){
        printf("search for source %u is not correct\n", s);
        return false;
      }
    }
    printf("PASS! multi-search verify\n");
    return true;
  }
};