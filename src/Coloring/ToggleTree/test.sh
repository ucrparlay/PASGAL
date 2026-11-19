#!/usr/bin/env bash
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
source ../../toggle_tree_config.sh
RUN_TEST_GRAPH "${TOGGLE_TREE_ROOT}" true void _sym.bin
