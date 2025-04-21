# PASGAL
PASGAL: Parallel And Scalable Graph Algorithm Library.  

## Prerequisite
+ g++ or clang with C++17 features support (tested with g++ 13.2.1 and clang 14.0.6) on Linux machines.  
+ We use [ParlayLib](https://github.com/cmuparlay/parlaylib) to support fork-join parallelism and some parallel primitives. It is provided as a submodule in our repository.  

## Algorithms
We include the following four algorithms in our repository. The source code can be found under ``src/``.  
* BFS: Breadth-First Search.  
* SSSP: Single-Source Shortest Paths. The SSSP algorithms are from the paper [[1]](#1).  
* BCC: Biconnected Components. The BCC algorithm is from the paper [[2]](#2).  
* SCC: Strongly Connected Components. The SCC algorithm is from the paper [[3]](#3).  
* basic_analytics: For computing no of vertices, edges, min_degree, max_degree, zero_degree_count.

## Compilation
A Makefile is provided in each subdirectory. For example, to compile BFS:  
```bash
cd src/BFS  
make  
```

## Running the code
Instructions on running the code will be provided when running the executables without any command line options. A sample output from BFS:  
> Usage: ./bfs [-i input_file] [-s] [-v]  
> Options:  
>         -i,     input file path  
>         -s,     symmetrized input graph  
>         -v,     verify result  

Graph Formats
--------
The graphs tested in our papers can be found at: https://pasgal-bs.cs.ucr.edu/.
The application can auto-detect the format of the input graph based on the suffix of the filename. Here is a list of supported graph formats: 
+ `.bin` The binary graph format from [GBBS](https://github.com/ParAlg/gbbs). It uses the compressed sparse row (CSR) format and is organized as follows:  
    + $n$ - number of vertices (64-bit variable)  
    + $m$ - number of edges (64-bit variable)  
    + sizes - sizes of this file in bytes, which is equal to $3\times8+(n+1)\times8+m\times4$ (64-bit variable)  
    + offset[] - offset[ $i$ ] (inclusive) and offset[ $i+1$ ] (exclusive) represents the range of neighbors list of the $i$-th vertex in the edges array (64-bit array of length $n+1$)  
    + edges[] - edges list (32-bit array of length $m$)  
+ `.adj` The adjacency graph format from [Problem Based Benchmark suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html).  

## Running Examples  

The graphs used in the [paper](#references) are available [here](https://pasgal-bs.cs.ucr.edu/bin/). They are in binary format, with `_sym` indicating undirected graphs, while others are directed.  

For demonstration, we will use:  
- **Directed graph**: `soc-LiveJournal1.bin`  
- **Undirected graph**: `soc-LiveJournal1_sym.bin`  

Since these graphs are large, a stack overflow may occur. To prevent this, run:  
```sh
ulimit -s unlimited
```

#### Running BFS
After compilation, three executables will be available in `src/BFS`:

* `bfs`: Runs our BFS algorithm from 5 random sources, repeating 6 times (ignoring the first warmup round).
* `seq-bfs`: Runs the standard sequential BFS in the same manner. A specific source can be set using `-r source`.
* `bfs_test`: Runs both our parallel BFS and the sequential BFS from 10 random sources, recording the average runtime in `bfs.tsv` and `seq-bfs.tsv`.


Run on an **undirected** graph:
```sh
./bfs -i path_to_graph/soc-LiveJournal1_sym.bin -s 
```
Run on a **directed** graph
```sh
./bfs -i path_to_graph/soc-LiveJournal1.bin
```

#### Running BCC
After compilation, three executables will be available in  `src/BCC`:
* `fast-bcc`:  Our parallel BCC algorithm [[2]](#2)
* `hopcroft-tarjan`:  Standard sequential BCC algorithm [[HT73]](https://dl.acm.org/doi/10.1145/362248.362272).
* `tarjan-vishkin`: A theoretically-efficient parallel BCC baseline [[TV85]](https://doi.org/10.1137/0214061).

BCC is only defined for **undirected** graphs:
```sh
./fast-bcc -i path_to_graph/soc-LiveJournal1_sym.bin -s 
```

#### Running SCC
After compilation, two executables will be avaible in  `src/SCC`:
* `scc`: Our parallel SCC algorithm [[3]](#3)
* `tarjan`: Standard sequential SCC algorithm. 

SCC is only defined for **directed** graphs.
```sh
./scc -i path_to_graph/soc-LiveJournal1.bin
```

#### Running SSSP
The SSSP algorithm supports **weighted** graphs, which are stored in the adjacency graph format and available [here](https://pasgal-bs.cs.ucr.edu/pbbs/)

After compilation, two executables will be avaible in `src/SSSP`:
* `sssp`: Our parallel SSSP algorithm [[1]](#1)
* `dijkstra`: Standard sequential SSSP algorithm. 

```sh
./sssp -i path_to_graph/soc-LiveJournal1_wgh18.adj
```

#### Running basic_analytics
After compilation, single executable will be available in `src/basic_analytics`. Run the executable as shown below:
```sh
./basic_analytics -i path_to_graph/soc-LiveJournal1_wgh18.adj
```


## References

Xiaojun Dong, Yan Gu, Yihan Sun, Letong Wang. [Brief Announcement: PASGAL: Parallel And Scalable Graph Algorithm Library](https://dl.acm.org/doi/10.1145/3626183.3660258). In *ACM Symposium on Parallelism in Algorithms and Architectures (SPAA)*, pp. 439-442, 2024.  

Xiaojun Dong, Yan Gu, Yihan Sun, Letong Wang. [PASGAL: Parallel And Scalable Graph Algorithm Library](https://arxiv.org/abs/2404.17101). *arXiv preprint: 2404.17101*, 2024.  

If you use our code, please cite our paper:
```
@inproceedings{dong2024pasgal,
  title = {Brief Announcement: PASGAL: Parallel And Scalable Graph Algorithm Library},
  author = {Dong, Xiaojun and Gu, Yan and Sun, Yihan and Wang, Letong},
  booktitle = {ACM Symposium on Parallelism in Algorithms and Architectures (SPAA)},
  year = {2024},
}
```

<a id="1">[1]</a>
Xiaojun Dong, Yan Gu, Yihan Sun, and Yunming Zhang. 2021. Efficient Stepping Algorithms and Implementations for Parallel Shortest Paths. In ACM Symposium on Parallelism in Algorithms and Architectures (SPAA). 184–197.

<a id="2">[2]</a> 
Xiaojun Dong, Letong Wang, Yan Gu, and Yihan Sun. 2023. Provably Fast and Space-Efficient Parallel Biconnectivity. In ACM Symposium on Principles and Practice of Parallel Programming (PPOPP). 52–65.

<a id="3">[3]</a>
Letong Wang, Xiaojun Dong, Yan Gu, and Yihan Sun. 2023. Parallel Strong Connectivity Based on Faster Reachability. ACM SIGMOD International Conference on Management of Data (SIGMOD) 1, 2 (2023), 1–29.


