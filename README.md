# PASGAL
PASGAL: Parallel And Scalable Graph Algorithm Library.  

## Repository Layout
PASGAL collects **multiple parallel implementations of the same graph problems**. Each problem has its own directory under `src/`, and each subdirectory inside it is one self-contained implementation following a different design approach:

```
src/
├── <Problem>/
│   ├── <Implementation A>/   # self-contained: own build & run scripts
│   └── <Implementation B>/
├── utils/                    # graph generators / converters / checkers
├── basic_analytics/          # basic graph statistics
└── graph.h, utils.h, ...     # headers shared by the PASGAL implementations
```

Implementations are developed and maintained independently, so each one keeps its own build system, command-line interface and output format. They are peers: no implementation is a baseline for another.

### Problems × Implementations

| Problem | [PASGAL](#pasgal-implementations) | [ToggleTree](#toggletree-implementations) |
|---|---|---|
| BFS: Breadth-First Search | `src/BFS/PASGAL` | `src/BFS/ToggleTree-dist` (distances)<br>`src/BFS/ToggleTree-parent` (BFS tree) |
| SSSP: Single-Source Shortest Paths | `src/SSSP/PASGAL` [[1]](#1) | `src/SSSP/ToggleTree-BellmanFord`<br>`src/SSSP/ToggleTree-Dijkstra` |
| BCC: Biconnected Components | `src/BCC/PASGAL` [[2]](#2) | |
| SCC: Strongly Connected Components | `src/SCC/PASGAL` [[3]](#3) | |
| KCore: k-Core Decomposition | `src/KCore/PASGAL` | `src/KCore/ToggleTree` |
| Coloring: Graph Coloring | | `src/Coloring/ToggleTree` (largest-first) |

An implementation may provide several variants of one problem; they are named `<Implementation>-<Variant>`.

### Implementations
* **PASGAL**: the implementations from the PASGAL papers ([References](#references)), built on [ParlayLib](https://github.com/cmuparlay/parlaylib) for fork-join parallelism and parallel primitives. Each directory also contains the standard sequential algorithm(s) for comparison.
* **ToggleTree**: implementations built on [Toggle Tree](https://github.com/ucrparlay/toggle-tree), a header-only parallel algorithm library based on the Toggle Tree data structure, with its own lightweight scheduler.

## Prerequisite
+ g++ or clang with C++20 support on Linux machines.  
+ Libraries are provided as git submodules under `external/` (run `git submodule update --init`):
    + [ParlayLib](https://github.com/cmuparlay/parlaylib), used by the PASGAL implementations.
    + [Toggle Tree](https://github.com/ucrparlay/toggle-tree), used by the ToggleTree implementations.

## PASGAL Implementations
A Makefile is provided in each implementation directory. For example, to compile the PASGAL BFS:  
```bash
cd src/BFS/PASGAL  
make  
```

### Running the code
Instructions on running the code will be provided when running the executables without any command line options. A sample output from BFS:  
> Usage: ./bfs [-i input_file] [-s] [-v]  
> Options:  
>         -i,     input file path  
>         -s,     symmetrized input graph  
>         -v,     verify result  

### Graph Formats
The graphs tested in our papers can be found at: https://pasgal-bs.cs.ucr.edu/.
The application can auto-detect the format of the input graph based on the suffix of the filename. Here is a list of supported graph formats: 
+ `.bin` The binary graph format from [GBBS](https://github.com/ParAlg/gbbs). It uses the compressed sparse row (CSR) format and is organized as follows:  
    + $n$ - number of vertices (64-bit variable)  
    + $m$ - number of edges (64-bit variable)  
    + sizes - sizes of this file in bytes, which is equal to $3\times8+(n+1)\times8+m\times4$ (64-bit variable)  
    + offset[] - offset[ $i$ ] (inclusive) and offset[ $i+1$ ] (exclusive) represents the range of neighbors list of the $i$-th vertex in the edges array (64-bit array of length $n+1$)  
    + edges[] - edges list (32-bit array of length $m$)  
+ `.adj` The adjacency graph format from [Problem Based Benchmark suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html).  

### Running Examples  

The graphs used in the [paper](#references) are available [here](https://pasgal-bs.cs.ucr.edu/bin/). They are in binary format, with `_sym` indicating undirected graphs, while others are directed.  

For demonstration, we will use:  
- **Directed graph**: `soc-LiveJournal1.bin`  
- **Undirected graph**: `soc-LiveJournal1_sym.bin`  

Since these graphs are large, a stack overflow may occur. To prevent this, run:  
```sh
ulimit -s unlimited
```

#### Running BFS
After compilation, three executables will be available in `src/BFS/PASGAL`:

* `bfs`: Runs the PASGAL BFS algorithm from 5 random sources, repeating 6 times (ignoring the first warmup round).
* `seq-bfs`: Runs the standard sequential BFS in the same manner. A specific source can be set using `-r source`.
* `bfs_test`: Runs both the PASGAL parallel BFS and the sequential BFS from 10 random sources, recording the average runtime in `bfs.tsv` and `seq-bfs.tsv`.


Run on an **undirected** graph:
```sh
./bfs -i path_to_graph/soc-LiveJournal1_sym.bin -s 
```
Run on a **directed** graph
```sh
./bfs -i path_to_graph/soc-LiveJournal1.bin
```

#### Running BCC
After compilation, three executables will be available in `src/BCC/PASGAL`:
* `fast-bcc`:  The PASGAL parallel BCC algorithm [[2]](#2)
* `hopcroft-tarjan`:  Standard sequential BCC algorithm [[HT73]](https://dl.acm.org/doi/10.1145/362248.362272).
* `tarjan-vishkin`: A theoretically-efficient parallel BCC baseline [[TV85]](https://doi.org/10.1137/0214061).

BCC is only defined for **undirected** graphs:
```sh
./fast-bcc -i path_to_graph/soc-LiveJournal1_sym.bin -s 
```

#### Running SCC
After compilation, two executables will be available in `src/SCC/PASGAL`:
* `scc`: The PASGAL parallel SCC algorithm [[3]](#3)
* `tarjan`: Standard sequential SCC algorithm. 

SCC is only defined for **directed** graphs.
```sh
./scc -i path_to_graph/soc-LiveJournal1.bin
```

#### Running SSSP
The SSSP algorithm supports **weighted** graphs, which are stored in the adjacency graph format and available [here](https://pasgal-bs.cs.ucr.edu/pbbs/)

After compilation, two executables will be available in `src/SSSP/PASGAL`:
* `sssp`: The PASGAL parallel SSSP algorithm [[1]](#1)
* `dijkstra`: Standard sequential SSSP algorithm. 

```sh
./sssp -i path_to_graph/soc-LiveJournal1_wgh18.adj
```

#### Running KCore
After compilation, a single executable `kcore` will be available in `src/KCore/PASGAL`. Run it without options to see its usage.

## Tools
#### basic_analytics
Computes the number of vertices and edges, and the min/max/zero in- and out-degree counts. Compile with `make` in `src/basic_analytics`, then run:
```sh
./basic_analytics -i path_to_graph/soc-LiveJournal1_wgh18.adj
```

#### utils
`src/utils` contains graph generators (`generate_random_graph`, `generate_grid_graph`), format converters (`pbbs2bin`, `symmetrize`), a graph checker (`validator`) and a diameter estimator (`get_diameter`). Compile with `make` in `src/utils` (`validator` is built with `make validator`).

## ToggleTree Implementations
Each ToggleTree implementation directory contains a `test.cpp` driver and two scripts:
* `test.sh`: compiles `test.cpp` and runs it on a single graph.
* `bench.sh`: compiles `test.cpp` and runs it on every graph in the benchmark list.

Both scripts can be launched from any working directory, and write `benchmark.csv` (running time) and `verify.csv` (result hash for a quick sanity check) next to the driver.

### Configuration
All ToggleTree implementations share one configuration file, `src/toggle_tree_config.sh`. It loads the runner functions provided by the Toggle Tree library, then sets the graph directory (`GRAPH_DIR`), the number of rounds, the compiler flags, the graph used by `test.sh` (`TEST`) and the graph list used by `bench.sh` (`GRAPHS`).

For machine-specific settings, create `src/toggle_tree_config.local.sh` (ignored by git); it is loaded last and overrides the shared values, e.g.:
```sh
GRAPH_DIR=/path/to/your/graphs
```

### Running
Every ToggleTree implementation is launched the same way, e.g.:
```sh
bash src/BFS/ToggleTree-dist/test.sh
bash src/SSSP/ToggleTree-Dijkstra/bench.sh
```
BFS, KCore and Coloring run on unweighted undirected graphs (`<graph>_sym.bin`); SSSP runs on weighted undirected graphs (`<graph>_sym_wgh16.bin`).

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


