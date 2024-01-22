#pragma once
#include <climits>

#include "graph.h"
#include "hashbag.h"
#include "reach.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "utils.h"

#include "multi_reach.h"
using namespace std;
using namespace parlay;

#define LOCAL

template <class Graph>
class SCC {
  using NodeId = typename Graph::NodeId;
  using EdgeId = typename Graph::EdgeId;
  using LabelT = size_t;

  static constexpr NodeId UINT_N_MAX = std::numeric_limits<NodeId>::max();

  // If top_bit is set, the vertex is not settle down
  static constexpr size_t TOP_BIT = size_t(UINT_N_MAX) + 1;
  static constexpr size_t VAL_MASK = UINT_N_MAX;
  static constexpr std::pair<NodeId, NodeId> EMPTY=std::pair(UINT_N_MAX,UINT_N_MAX);

  const Graph &G;
  const Graph &GT;

  LabelT label_offset;

 public:
  SCC() = delete;
  SCC(const Graph& _G, const Graph& _GT) : G(_G), GT(_GT){}

  size_t trim1(sequence<LabelT>& label){
    // the vertices whose in-degree or out-degree is zero
    auto zeros = parlay::pack_index(
      parlay::tabulate(G.n, [&](NodeId i){
      return (G.offsets[i+1]-G.offsets[i]==0) || (GT.offsets[i+1]-GT.offsets[i]==0);}));
    parallel_for(0, zeros.size(),[&](NodeId i){
      NodeId v = zeros[i];
      label[v]=(size_t) i|TOP_BIT;
    });
    return zeros.size();
  }
  size_t scc_first(sequence<LabelT>&label){
    // auto degree_product = parlay::tabulate(G.n, [&](NodeId i){
    //   return size_t(G.offsets[i+1]-G.offsets[i])*\
    //           size_t(GT.offsets[i+1]-GT.offsets[i]);
    // });
    auto identities = parlay::delayed_tabulate(G.n,[&](NodeId i){return i;});
    auto NON_ZEROS=parlay::filter(identities, [&](NodeId v){return !(label[v]&TOP_BIT);});
    auto P = parlay::random_shuffle(NON_ZEROS);
    NodeId source = P[0];
    // auto max_candidate = parlay::max_element(degree_product);
    // NodeId source = max_candidate - degree_product.begin();
    printf("scc_first source: %u\n", source);
    REACH<Graph> reach_solver(G.n);
    // TODO: modify single reach to skip ZEROS
    parlay::sequence<bool> forward_reach(G.n);
    parlay::sequence<bool> backward_reach(G.n);
    reach_solver.reach(source, G, GT, forward_reach);
    reach_solver.reach(source, GT, G, backward_reach);
    parallel_for(0, G.n, [&](NodeId i){
      if (forward_reach[i]&& backward_reach[i]){
        label[i]=label_offset | TOP_BIT;
      }else if (!(TOP_BIT&label[i]) && (forward_reach[i]||backward_reach[i])){
        label[i]=label_offset;
      }
    });
    return parlay::reduce(parlay::tabulate<size_t>(G.n, [&](NodeId i){
      return (size_t)(forward_reach[i]&& backward_reach[i]);}));
  }
  void scc(sequence<LabelT>& label){
    parlay::internal::timer t("scc");
    parlay::parallel_for(0, G.n,[&](size_t i){label[i]=0;});
    label_offset=0;
    size_t n_zeros=trim1(label);
    t.next("trim1");
    printf("n_zeros: %ld\n", n_zeros);
    label_offset+=n_zeros;
    size_t n_scc1 = scc_first(label);
    label_offset+=1;
    printf("n_scc1: %ld\n", n_scc1);
    if (n_scc1 < G.n/100000){
      n_scc1 = scc_first(label);
      printf("n_scc2: %ld\n",n_scc1);
      label_offset+=1;
    }
    t.next("first round");
    auto P = parlay::random_permutation((NodeId)G.n);
    auto vertices = parlay::filter(P, [&](NodeId i){return !(TOP_BIT&label[i]);});
    // printf("remaining vertices: %ld\n", vertices.size());
    t.next("permute & filter");
    // multi_reach
    float beta = 1.5;
    size_t step = 2;
    size_t start = 0;
    size_t end=start+step;
    size_t round = 0;
    // MULTI_REACH<Graph,LabelT> multi_reach_solver(vertices.size(), G.n);
    MULTI_REACH<Graph,LabelT> multi_reach_solver(G.n, G.n);
    size_t fwd_m=1; size_t bwd_m=1;
    size_t n_remain=vertices.size();
    printf("n_remain %zu\n", n_remain);
    // size_t kCacheLineSz=128;
    // parlay::sequence<size_t> cts(parlay::num_workers()*kCacheLineSz);
    // for (size_t i = 0;i<parlay::num_workers(); i++){cts[i*kCacheLineSz]=0;}
    // parlay::internal::timer t1("  ",false);
    // float multi_search_time=0;
    // float set_intersect_time=0;
    parlay::internal::timer t_round(" ",false);
    while (end < n_remain){
      round++;
      end = std::min(start+step, n_remain);
      // printf("Round: %ld  start: %ld  end: %ld n_remain: %lu\n", round, start, end, n_remain);
      auto sources = parlay::filter(vertices.cut(start, end), 
            [&](NodeId i){return !(label[i]&TOP_BIT);});
      // fwd_m=std::max((size_t)floor(beta*fwd_m),2*n_remain);
      // bwd_m=std::max((size_t)floor(beta*bwd_m),2*n_remain);
      fwd_m=2*max((size_t)min((NodeId)ceil(0.3*n_remain),(NodeId)6000000),
                (size_t)(beta)*fwd_m);
      bwd_m=2*max((size_t)min((NodeId)ceil(0.3*n_remain),(NodeId)6000000),
                (size_t)(beta)*bwd_m);
      // fwd_m=2*max((size_t)ceil(0.3*n_remain),(size_t)(beta)*fwd_m);
      // bwd_m=2*max((size_t)ceil(0.3*n_remain),(size_t)(beta)*bwd_m);
      // printf("  frontier: %ld fwd_m: %ld bwd_m: %ld\n", sources.size(), fwd_m, bwd_m);
      hash_table<NodeId, NodeId> fwd_table(fwd_m, EMPTY);
      hash_table<NodeId, NodeId> bwd_table(bwd_m, EMPTY);
      // t1.start();
      multi_reach_solver.multi_reach_safe(G,label,sources,fwd_table,label_offset);
      multi_reach_solver.multi_reach_safe(GT,label,sources,bwd_table,label_offset);
      // multi_search_time+=t1.stop();

      fwd_m=fwd_table.size();
      bwd_m=bwd_table.size();
      // printf("forward_m %ld backward_m %ld\n", fwd_m, bwd_m);

      // t1.start();
      auto& smaller_t = (fwd_m <= bwd_m) ? fwd_table : bwd_table;
      auto& larger_t = (fwd_m > bwd_m) ? fwd_table : bwd_table;

      auto intersection_f = [&](auto& kv) {
        NodeId key = std::get<0>(kv);
        LabelT value = std::get<1>(kv);
        // if (larger_t.contain(kv)) {
        //   // write_max(&labels[key],value);
        //   bool suc = 0; LabelT c;
        //   do c = *(&label[key]);
        //   while ((value>c) && !(suc = atomic_compare_and_swap(&label[key], c, value)));
        //   if (suc && ((c&TOP_BIT)>0)){cts[parlay::worker_id()*kCacheLineSz]++;}
        // } else {
        //   write_max(&label[key],value|TOP_BIT);
        // }
        if (larger_t.contain(kv)){
          write_max(&label[key], value|TOP_BIT);
        }else{
          write_max(&label[key], value);
        }
      };
      smaller_t.map(intersection_f);
      // size_t n_processed = 0;
      // for(size_t i = 0; i<parlay::num_workers(); i++){
      //   n_processed+=cts[i*kCacheLineSz];
      //   cts[i*kCacheLineSz]=0;
      // }
      // printf("  fwd_m: %ld bwd_m: %ld\n",fwd_m, bwd_m);
      auto set_label=[&](auto& kv){
        NodeId key=std::get<0>(kv);
        LabelT value = std::get<1>(kv);
        write_max(&label[key], value);
      };
      larger_t.map(set_label);
      // set_intersect_time += t1.stop();
      // printf("  end larger_map for setting labels\n");
      label_offset += sources.size();
      start=end;
      step = floor(beta*step);
      t_round.next("  round");
    }
    // printf("  multi-search time: %f\n", multi_search_time);
    // printf("  set-intersection time: %f\n", set_intersect_time);
  }
  void status(sequence<LabelT>& label){
    auto sorted_label = parlay::sort(label);
    auto flag = parlay::tabulate(G.n,[&](size_t i){
      if (i==0 || sorted_label[i]!= sorted_label[i-1]){
        return true;
      }else{
        return false;
      }
    });
    auto index = parlay::pack_index(flag);
    index.push_back(G.n);
    auto max_scc = parlay::max_element(parlay::delayed_tabulate(index.size()-1,[&](size_t i){return index[i+1]-index[i];}));
    printf("number of scc: %zu\n", index.size()-1);
    printf("max scc: %zu\n", *max_scc);
  }
};

