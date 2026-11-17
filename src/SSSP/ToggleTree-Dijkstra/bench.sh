#!/usr/bin/env bash
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
source ../../toggle_tree_config.sh
RUN_BENCH_GRAPH "${TOGGLE_TREE_ROOT}" true uint32_t _sym_wgh16.bin
