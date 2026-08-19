#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

Rscript -e 'files <- list.files("replication/code/r", pattern="[.]R$", full.names=TRUE); for (f in files) parse(f); cat("PASS: all R files parse.\n")'
Rscript replication/code/r/smoke_test_class_sampling.R

OUTPUT_DIR="$(mktemp -d "${TMPDIR:-/tmp}/representation-multiplicity-equality.XXXXXX")"
trap 'rm -rf "$OUTPUT_DIR"' EXIT
RI_OUTPUT_DIR="$OUTPUT_DIR" Rscript replication/code/r/verify_reported_equality.R
python3 - "$OUTPUT_DIR/class_sampled_reported_equality.csv" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    row = next(csv.DictReader(handle))
for field in (
    "maximum_absolute_prediction_difference",
    "maximum_absolute_variance_difference",
    "maximum_absolute_standard_error_difference",
):
    if float(row[field]) != 0:
        raise SystemExit(f"FAIL: {field} is nonzero")
print("PASS: finite-array predictions, variances, and standard errors are identical.")
PY
