#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

Rscript replication/code/r/run_dimension_scaling.R
Rscript replication/code/r/run_sign_disagreement.R
Rscript replication/code/r/run_seed_sensitivity.R
RI_RUN_SCOPE=simulation Rscript replication/code/r/run_class_sampling_study.R

printf 'Simulation reruns written under %s/rerun_outputs/.\n' "$ROOT"
