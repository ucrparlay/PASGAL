# Reachability

Parallel reachability (transitive closure of a single source) on the PASGAL
graph format.  The implementation is a **sparse-only push BFS** with a
multi-hop local queue.  Direction optimization (sparse↔dense switching) is
intentionally out of scope — see the plan section below.

## Files

| File | Purpose |
|---|---|
| `reachability.h` | Algorithm: `Reachability<Graph>` class (header-only template). Hashbag frontier. |
| `reachability.cpp` | CLI driver, timing harness, sequential-BFS verifier. |
| `reachability_winning_tree.{h,cpp}` | Same algorithm with `WinningTree` instead of hashbag. Kept for reference (slower; see results below). |
| `run_reachability.py` | Single-config sweep over the standard graph set. |
| `sweep_reachability.py` | Parameter sweep, **resumable** (skips configs already in the output TSV). |
| `Makefile` / `CMakeLists.txt` | Build. Builds both `reachability` and `reachability_winning_tree`. |

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
| `beta` (`-b`) | Edges per burst (also the `parallel_for` granularity in `visit_neighbors_parallel` and the `deg < beta` threshold for sequential vs parallel fan-out). | **2048** |
| `max_queue_size` (`-q`) | Vertex cap on the per-task local queue. | **2000** |
| `mode` (`-t`) | Burst exit condition: `dual` / `vertex` / `edge`. | **`vertex`** |

Inside each burst:

```cpp
const size_t queue_size = min(max_queue_size,
                              max(1, num_threads * beta / frontier_size));
const size_t max_edges  = queue_size * EDGE_PER_CACHELINE;

while (front < rear && /* mode-dependent condition */) { ... }
//   mode = vertex : vertices_visited < max_queue_size
//   mode = edge   : edges_processed   < max_edges
//   mode = dual   : both
```

`queue_size` is the **adaptive** per-task cap — small when the frontier is
large (so there's enough work for everyone), large when the frontier is
small (so each task chain-walks further).

The local queue is a `thread_local std::vector<NodeId>`, lazily resized to
`max_queue_size` on first use, reused across tasks on the same worker.

### Implementation notes

- **`visited`** is `parlay::sequence<std::atomic<bool>>` with
  `exchange(true, std::memory_order_relaxed)` on the hot path.  Relaxed
  ordering is sufficient because `visited[v]` is just a dedup flag — no
  happens-before relation needs to be established between vertices.
- **No `in_frontier` array.**  Every caller of `add_to_frontier(v)` has
  just won `visited[v].exchange(true)`, so each `v` is uniquely owned by
  one thread per round; the bag insert needs no extra dedup.
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
               [-b <beta>] [-q <max_queue_size>] [-t <mode>]
```

| Flag | Meaning |
|---|---|
| `-i` | Input graph path (binary CSR, see `graph.h`). Required. |
| `-o` | Output TSV path (default `reachability.tsv` in cwd). |
| `-s` | Treat graph as symmetric (skip `make_inverse`). |
| `-v` | Verify each result against sequential BFS. |
| `-r <id>` | Single source; otherwise averages over 5 hash-derived sources. |
| `-b <n>` | Burst edge budget (default 2048). |
| `-q <n>` | Local-queue vertex cap (default 2000). |
| `-t <m>` | Threshold mode `dual` / `vertex` / `edge` (default `vertex`). |

Output TSV columns: `graph\tsource\ttime\tbeta\tmax_queue_size\tmode`.

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

## Empirical results

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

## Plan (sparse-only refinement)

The current plan, kept at `~/.claude/plans/okay-good-do-you-cached-porcupine.md`,
is sparse-only.  Direction optimization, AST-based connectivity for
symmetric graphs, multi-source bitmask, and continuation-style spawning
are all explicitly out of scope under the current paper framing.

| Step | Status | Notes |
|---|---|---|
| 1. Drop `in_frontier`, keep atomic `visited`, relax memory ordering | **Done** | One fewer per-node atomic; ~5–10% speedup expected on small-frontier workloads. |
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
