#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${1:-$ROOT/replication/data/kenya}"

cd "$ROOT"
RI_SOURCE_DIR="$DATA_DIR" Rscript replication/code/r/run_fixed_budget_study.R
RI_SOURCE_DIR="$DATA_DIR" RI_RUN_SCOPE=all Rscript replication/code/r/run_class_sampling_study.R

printf 'Kenya reruns written under %s/rerun_outputs/.\n' "$ROOT"
