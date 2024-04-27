# PASGAL
PASGAL: Parallel And Scalable Graph Algorithm Library.  

Prerequisite
--------
+ g++ or clang with C++17 features support (tested with g++ 13.2.1 and clang 14.0.6) on Linux machines.  
+ We use [ParlayLib](https://github.com/cmuparlay/parlaylib) to support fork-join parallelism and some parallel primitives. It is provided as a submodule in our repository.  

Algorithms
--------
We include the following four algorithms in our repository. The source code can be found under ``src/``.  
* BFS: Breadth-First Search  
* BCC: Biconnected Components.  The BCC algorithm is from the paper [[1]](#1).
* SCC: Strongly Connected Components.  The SCC algorithm is from the paper [[2]](#2).
* SSSP: Single-Source Shortest Paths  


Compilation
--------
A Makefile is provided in each subdirectory. For example, to compile BFS:  
```bash
cd src/BFS  
make  
```

Running the code
--------
Instructions on running the code will be provided when running the executables without any command line options. A sample output from BFS:  
> Usage: ./bfs [-i input_file] [-s] [-v]  
> Options:  
>         -i,     input file path  
>         -s,     symmetrized input graph  
>         -v,     verify result  

Graph Formats
--------
The application can auto-detect the format of the input graph based on the suffix of the filename. Here is a list of supported graph formats: 
+ `.bin` The binary graph format from [GBBS](https://github.com/ParAlg/gbbs). It uses the compressed sparse row (CSR) format and is organized as follows:  
    + $n$ - number of vertices (64-bit variable)  
    + $m$ - number of edges (64-bit variable)  
    + sizes - sizes of this file in bytes, which is equal to $3\times8+(n+1)\times8+m\times4$ (64-bit variable)  
    + offset[] - offset[ $i$ ] (inclusive) and offset[ $i+1$ ] (exclusive) represents the range of neighbors list of the $i$-th vertex in the edges array (64-bit array of length $n+1$)  
    + edges[] - edges list (32-bit array of length $m$)  
+ `.adj` The adjacency graph format from [Problem Based Benchmark suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html).  

## References
<a id="1">[1]</a> 
Dong, Xiaojun, et al. "Provably Fast and Space-Efficient Parallel Biconnectivity." Proceedings of the 28th ACM SIGPLAN Annual Symposium on Principles and Practice of Parallel Programming. 2023.

<a id="2">[2]</a>
Wang, Letong, et al. "Parallel Strong Connectivity Based on Faster Reachability." Proceedings of the ACM on Management of Data 1.2 (2023): 1-29.
