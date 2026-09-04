#!/bin/sh
set -eu
cd "$(dirname "$0")/../paper/source"
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
