# Reachability

Parallel reachability (transitive closure of a single source) on the PASGAL
graph format.  Push (sparse) relaxation with a multi-hop local queue, pull
(dense) relaxation over a frontier bitmap, and direction optimization
switching between them by frontier density.  The frontier is a chunk bag.

## Files

| File | Purpose |
|---|---|
| `reachability.h` | Algorithm: `Reachability<Graph>` class (header-only template). Chunkbag frontier, push/pull direction optimization. |
| `reachability.cpp` | CLI driver, timing harness, sequential-BFS verifier. |
| `reachability_bench.cpp` | Beta sweep across many sources per invocation; reports time, rounds, dense rounds. |
| `reachability_winning_tree.{h,cpp}` | The *previous* sparse-only algorithm with `WinningTree` instead of a bag. Kept for reference (slower; see results below). |
| `Makefile` / `CMakeLists.txt` | Build. `make` builds `reachability`, `reachability_bench`, `reachability_winning_tree`. |

The driver script `run_reachability.py` lives in `Synchronization-Model/`,
not here.  The previous hashbag implementation and its sweep outputs are
kept outside the repo, in `~/Projects/Model/attic/reachability/`.

## Algorithm

Standard sparse-frontier BFS with a *local-queue* extension that lets each
parallel-for task chain-walk many BFS levels before returning to the
scheduler.

```
visited[s] = true; bag.insert(s)
while bag not empty:
  frontier <- bag.pack_into                 # snapshot of next frontier
  parallel_for f in frontier:
    local_queue = [f]
    while burst budget left:
      u = local_queue.pop()
      for v in out-neighbors(u):
        if visited[v].exchange(true) == false:
          if local_queue.has_space():
            local_queue.push(v)
          else:
            bag.insert(v)
    bag.insert(remaining local_queue entries)
```

The burst budget caps each task's work to amortize fork/sync cost while
keeping the work-stealing scheduler effective.

### Local-queue parameters

| Parameter | Meaning | Default |
|---|---|---|
| `beta` (`-b`) | Edges per burst, a runtime member. | **2048** |
| `MAX_QUEUE` | Vertex cap on the per-task local queue, `static constexpr`. | **4096** |
| `SPARSE_TH` | Go dense when the frontier holds ≥ `n / SPARSE_TH` vertices. | **20** |
| `BLOCK_SIZE` | `deg <` this expands a vertex sequentially, else in parallel. | **1024** |

Each burst walks until the edge budget runs out:

```cpp
while (front < rear && edges < beta) { ... }
```

A flat budget replaced the earlier adaptive `queue_size` (which scaled with
`num_threads * beta / frontier_size`) and the `dual`/`vertex`/`edge` mode
selector; both were measured to be noise here.  The local queue is a plain
stack array of `MAX_QUEUE` entries rather than a `thread_local` vector.

### Implementation notes

- **`visited`** is `parlay::sequence<uint8_t>` claimed with
  `compare_and_swap`.  It is only a dedup flag, so no happens-before
  relation needs to be established between vertices.
- **`in_frontier` is back.**  It was dropped when the algorithm was
  sparse-only, because every caller of `add_to_frontier(v)` had just won
  the `visited[v]` CAS and so uniquely owned `v`.  That no longer holds:
  `dense_to_sparse()` refills the bag from the bitmap, whose vertices had
  their `visited` CAS won earlier inside `dense_relax()`, so the bag insert
  needs its own dedup flag.
- **The pull pass exits on its first hit.**  With no distance to minimize,
  the first reachable in-neighbor settles a vertex, so the in-edge scan
  always breaks early — unlike delta-stepping BFS, whose scan must stay
  permissive.  For the same reason the local-queue walk costs no
  speculation, so `beta` trades rounds against redundant work only.
- **Verification** uses `BFS/seq-bfs.h` (a plain `std::queue` BFS) as
  ground truth — no parlay primitives, fully independent.

## Build

```bash
make                       # default flags: clang++ -O3 -march=native -pthread
make OPENCILK=1            # use OpenCilk
make SERIAL=1              # single-threaded debug build
```

## Run

### Single-config

```bash
./reachability -i <graph.bin> [-o <output.tsv>] [-s] [-v] [-r <source>] \
               [-b <beta>] [-D]
```

| Flag | Meaning |
|---|---|
| `-i` | Input graph path (binary CSR, see `graph.h`). Required. |
| `-o` | Output TSV path (default `reachability.tsv` in cwd). |
| `-s` | Treat graph as symmetric (skip `make_inverse`). |
| `-v` | Verify each result against sequential BFS. |
| `-r <id>` | Single source; otherwise averages over 5 hash-derived sources. |
| `-b <n>` | Burst edge budget (default 2048). |
| `-D` | Disable direction optimization (push only). |

Output TSV columns: `graph\tsource\ttime\tbeta\trounds\tdense_rounds`.

### Beta sweep

```bash
./reachability_bench -i <graph.bin> [-s] (-r <source>).. [-n <reps>] \
                     (-b <beta>).. [-D | -A]
```

Many sources per invocation, so a large graph is read and transposed once;
`-A` measures both direction arms in the same process.  Emits
`REACH_SOURCE` and `REACH_TIME <mode> <beta> <rep> <sec> <rounds>
<dense_rounds> <reached>`, and fails if any beta or arm reaches a different
vertex count than the first measurement for that source.

### Standard sweep (default config)

```bash
python3 run_reachability.py
```

One run per graph at the default parameters.  Output → `output/reachability.tsv`,
log → `logs/reachability_<date>.log`.

### Parameter sweep (post-cleanup)

```bash
python3 sweep_reachability.py
```

Sweeps `β ∈ {2048}` × `q ∈ {500, 1000, 2000, 5000, 10000}` × `mode ∈ {dual, vertex, edge}`
over the 12 standard graphs (5 sources each) — **180 timed runs + 180
sequential-BFS verifications**.  The script is resumable: existing rows
in the TSV are detected and skipped.

Output → `results/<date>_reachability_post_cleanup_sweep.tsv`.

## Empirical results (previous sparse-only implementation)

Everything in this section predates direction optimization and the chunkbag
frontier, and the `q`/`mode` knobs it tunes no longer exist.  Kept because
the per-graph rankings are still informative; the numbers are not current.

From the wide sweep at `results/2026-04-28_reachability_threshold_sweep.tsv`
(12 graphs × 3 modes × 16 q values × 5 sources):

- **`(mode=vertex, q=2000)` is the universal default.**
- Wide plateau: `q ∈ [1000, 20000]` is within 5% of best on geomean.
- Mode ranking: vertex sweeps the top-6 geomean slots; edge and dual tie
  ~3% behind.
- **Worst-case overhead vs per-graph oracle: 1.54×** (CHEM_5).
- Old default `(dual, 1000)` had worst-case 2.50× (Cosmo50).

Per-graph oracle vs default:

| graph | oracle (mode, q) | default vs oracle |
|---|---|---|
| africa.bin | (edge, 200) | 1.38× |
| north-america.bin | (vertex, 2000) | **1.00×** ✓ |
| asia.bin | (vertex, 200) | 1.27× |
| europe.bin | (vertex, 1000) | 1.01× |
| CHEM_5.bin | (vertex, 20000) | **1.54×** (worst) |
| GeoLifeNoScale_5.bin | (dual, 10) | 1.13× |
| GeoLifeNoScale_10.bin | (edge, 50000) | 1.01× |
| Cosmo50_5.bin | (vertex, 20000) | 1.25× |
| grid_1000_100000.bin | (vertex, 20000) | 1.24× |
| grid_1000_100000_03.bin | (edge, 2000) | 1.04× |
| hugetrace-00020_sym.bin | (vertex, 10000) | 1.24× |
| hugebubbles-00020_sym.bin | (vertex, 10000) | 1.25× |

Of the 12 oracles: **9 use vertex mode, 2 edge, 1 dual** — vertex is
consistently the right family.  Only the actual q varies (200–20000).

### Hashbag vs WinningTree

Tested separately: hashbag wins by ~3.5× geomean, 1.4–16× per graph.
Reason: WT's phase-concurrent insert-during-iterate machinery is overhead
without a level-correctness benefit (no delta-stepping in plain
reachability).  `reachability_winning_tree` exists for reproducibility
but is not the recommended path.

## Plan (superseded)

This plan was written when the algorithm was deliberately sparse-only.
Direction optimization has since been implemented, so the framing below no
longer holds; AST-based connectivity for symmetric graphs, multi-source
bitmask, and continuation-style spawning remain out of scope.

| Step | Status | Notes |
|---|---|---|
| 1. Drop `in_frontier`, keep atomic `visited`, relax memory ordering | **Reverted** | `dense_to_sparse()` reinserts vertices whose `visited` CAS was already won, so the bag needs its own dedup flag. |
| 2. Re-enable `-v` in `sweep_reachability.py` and confirm correctness post-cleanup | **Pending** | Script already updated; just needs to be run. |
| 3. Per-round diagnostic logging (paper instrumentation) | On hold | Round#, frontier size, sum of degrees, wall time — for α/β-validation plots. |
| 4. Theoretical write-up of `q ≈ 2000` from machine α/β | On hold | Calculation goes in the paper, not the code. |

### Reading-list grounding

Beyond the SCC paper (which the early plan leaned on), other prior work
contributes mostly framing rather than new algorithmic options:

- **PASGAL paper** (SPAA'24, [arxiv:2404.17101](https://arxiv.org/abs/2404.17101)):
  VGC = our local-queue multi-hop expansion; matches what we already have.
- **k-Core paper** (SIGMOD'25, [arxiv:2502.08042](https://arxiv.org/abs/2502.08042)):
  introduces "burdened span" cost model with ω ≈ 15 000 cycles for
  fork/join — direct precedent for the α/β framing.
- **SCC reachability** (SIGMOD'23, [arxiv:2303.04934](https://arxiv.org/abs/2303.04934)):
  canonical implementation of *direction-optimized* reachability.  Out of
  scope by user direction (sparse only).
- **Stepping/SSSP** (SPAA'21, [arxiv:2105.06145](https://arxiv.org/abs/2105.06145)):
  BFS is the δ=1 instance; LaB-PQ degenerates to a single-tier hashbag for
  unweighted reachability.  No new technique.
- **Biconnectivity** (PPoPP'23, [arxiv:2301.01356](https://arxiv.org/abs/2301.01356)):
  uses arbitrary spanning trees for connectivity (polylog span), not BFS.
  Different algorithm class; out of scope.
