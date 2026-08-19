#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
"$ROOT/paper/source/build.sh"

PDF="$ROOT/paper/representation_multiplicity_causal_forests.pdf"
if command -v pdfinfo >/dev/null 2>&1; then
  pages="$(pdfinfo "$PDF" | awk '/^Pages:/ {print $2}')"
  if [[ -z "$pages" || "$pages" -gt 45 ]]; then
    printf 'Paper page-count check failed: %s pages.\n' "${pages:-unknown}" >&2
    exit 1
  fi
  printf 'PASS: paper has %s pages (limit 45).\n' "$pages"
fi
