# Replication

## Check the Published Tables

From the repository root, run

```sh
python3 scripts/verify_replication.py
```

The verifier checks file hashes and reconstructs both tables from the included
CSV results. Python 3.10 or later is sufficient.

## Monte Carlo

Table 1 uses ten paired replications, 3,000 training observations, and a shared
independent evaluation array of 10,000 observations. Multiplicity is
5, 10, 40, or 400. Minimum node size is 5, 10, or 20. Each fit uses 2,000 trees
and `grf` 2.4.0. The three `mtry` rules are 2, the package default, and one third
of the number of columns rounded down. Class sampling applies these rules to
the six classes. The rows with multiplicity 1 contain no repeated columns.

The reference environment is R 4.6.1 with `grf` 2.4.0. Install that package and run

```sh
python3 scripts/run_monte_carlo.py --rscript /path/to/Rscript
```

The script runs one thread per forest and defaults to one worker per available
core minus one. Results go to `rerun_outputs/monte_carlo`. For a short run, use
`--replications 3 --multiplicities 5 --node-sizes 5 --trees 50 --n-train 200`.

Training data use seed `20261007 + r`, forest fits use seed
`20261007 + 30000 + r`, and the evaluation array uses seed `20261008`.
All inputs are rank-coded on the pooled training and evaluation array.

## Job Training and Earnings

Download `18515_Data_and_Programs.zip` from the
[Econometric Society](https://www.econometricsociety.org/publications/econometrica/2022/05/01/automatic-debiased-machine-learning-causal-and-structural).
Its SHA-256 is
`89b0e810319eec93bc40418ac802820d27c96d3c1bf6ea9786f6b75fffa169db`.
Extract the archive and locate `rrr_lasso_NSW_blackbox`.

The reference environment is Intel R 4.0.3 with the historical macOS binary
of `randomForest` 4.6-14 and `nnet` 7.3-14. Each outcome forest uses 1,000 trees,
bootstrap sampling, minimum node size 5, and five folds. Matching the R and
package binary builds is required to reproduce the original forest draws.

```sh
python3 scripts/run_empirical.py \
  --rscript /path/to/Intel-R-4.0.3/bin/Rscript \
  --package-root /path/to/rrr_lasso_NSW_blackbox
```

The script first runs the original estimator at seed 1. Paired comparisons
retain its samples, folds, Riesz estimates, and fold-level random states.
Seeds 2--5 change forest randomness while retaining the seed-1 folds.
Results go to `rerun_outputs/empirical`.

Specifications 1 and 2 remove columns with equal unordered threshold
partitions. Specification 3 removes eight known equivalent columns, reducing
344 columns to 336. These are the four degree-one polynomial terms and their
treatment interactions. Class sampling assigns these pairs identical ranks
before grouping. Retained raw values remain unchanged in the reduced-column
fits. The removed-column maps are included with the empirical results.

Table 2 reports ATET levels and standard errors. Brackets give signed changes
from the all-column fit. The seed benchmark reports mean absolute ATET
changes and mean SR, with maxima underneath. Eight original ATETs and all nine
standard errors match the published values to two decimal places. The CPS
specification 3 ATET is 1906.63 versus the published 1906.62.

## Build the Paper

```sh
./scripts/build_paper.sh
```

The self-contained source requires a LaTeX installation with `latexmk`.
