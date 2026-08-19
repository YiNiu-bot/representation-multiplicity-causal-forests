#!/usr/bin/env bash
set -euo pipefail

SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1787097600}"
export FORCE_SOURCE_DATE="${FORCE_SOURCE_DATE:-1}"
export TZ="${TZ:-UTC}"
cd "$SOURCE_DIR"
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
if grep -Eq 'Citation .* undefined|Reference .* undefined|There were undefined references|Overfull \\hbox' main.log; then
  printf 'LaTeX build contains unresolved citations, references, or overfull boxes.\n' >&2
  exit 1
fi
cp main.pdf ../representation_multiplicity_causal_forests.pdf
latexmk -c main.tex >/dev/null
rm -f main.pdf main.bbl
printf 'Built %s\n' "$SOURCE_DIR/../representation_multiplicity_causal_forests.pdf"
