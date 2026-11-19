#!/usr/bin/env bash
# PASGAL-side configuration for toggle-tree benchmarks.
#
# Usage (from a src/<ALGO>/ToggleTree/ script):
#   source "$(dirname "${BASH_SOURCE[0]}")/../../toggle_tree_config.sh"
#   RUN_TEST_GRAPH "${TOGGLE_TREE_ROOT}" true void _sym.bin
#
# The runner functions (COMPILE / RUN_TEST_GRAPH / RUN_BENCH_GRAPH / ...) come
# from upstream's example_of_config.sh; nothing inside external/ is modified.
# They read the variables below at call time, so overriding them here works.

PASGAL_SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOGGLE_TREE_ROOT="${PASGAL_SRC_DIR}/../external/toggle-tree"

# 1. Upstream runner functions + defaults
source "${TOGGLE_TREE_ROOT}/dataset/example_of_config.sh"

# 2. PASGAL settings (override upstream defaults)

# Your directory containing graphs
GRAPH_DIR=/data/graphs/bin

# Number of test rounds after warmup
NUM_ROUNDS=5

# Compiler and flags
CXX=g++
CXXFLAGS=(-O3 -std=c++20 -pthread -march=native -Wall)

# Uncomment to use parlaylib scheduler (path is relative to TOGGLE_TREE_ROOT).
# PARLAY_DIR=../parlaylib/include

# Clear page cache before each testcase.
CLEAR_CACHE=true

# The single graph run by test.sh
TEST=twitter

GRAPHS=(
    friendster
    twitter
    com-orkut
    socfb-konect
    soc-LiveJournal1
    # clueweb
    sd_arc
    arabic
    uk-2002
    indochina
    planet
    europe
    asia
    africa
    RoadUSA
    Cosmo50_5
    GeoLifeNoScale_10
    hugebubbles-00020
    hugetrace-00020
    CHEM_5
)

# 3. Optional per-machine overrides (not tracked by git)
if [[ -f "${PASGAL_SRC_DIR}/toggle_tree_config.local.sh" ]]; then
    source "${PASGAL_SRC_DIR}/toggle_tree_config.local.sh"
fi
