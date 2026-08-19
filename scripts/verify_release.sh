#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if command -v shasum >/dev/null 2>&1; then
  shasum -a 256 -c SHA256SUMS
elif command -v sha256sum >/dev/null 2>&1; then
  sha256sum -c SHA256SUMS
else
  printf 'No SHA-256 checksum command found.\n' >&2
  exit 1
fi

python3 scripts/validate_release.py
python3 replication/code/python/audit_tables.py

for check in replication/code/javascript/verify_*.js; do
  node "$check"
done

printf 'PASS: frozen release integrity, Tables I--III, and deterministic checks.\n'
