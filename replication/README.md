# Replication

## 1. Reconstruct the included tables

```sh
python3 scripts/verify_replication.py
```

Run from the repository root. Python 3.10+ is the only dependency. The
included CSVs contain all replication-level statistics used by the paper,
including 100 paired observations per Monte Carlo cell. The verifier checks
the manifest, computes means and Monte Carlo standard errors, reconstructs
paired loss differences, and checks all six native and canonical empirical
comparisons. It prints a JSON table reconstruction.

## 2. Refit the Monte Carlo

Reference environment: R 4.6.1, grf 2.4.0, Rcpp 1.1.2, Matrix 1.7-5,
15 threads per forest. Detailed session information is in `metadata/`.
Install grf 2.4.0 from the CRAN archive in an R library of your choice.
No R package or binary is redistributed here.

```sh
export RI_RLIB=/path/to/R-library
python3 scripts/reproduce.py --stage monte-carlo --rscript /path/to/Rscript
```

Defaults are 100 replications per dimension, 3,000 training observations,
30,000 independent evaluation observations, and 2,000 trees. The dimension
and seed-control experiments use separate recorded data-seed sequences.
Within each experiment, comparisons are paired. Class sampling is evaluated
at p=500 only. Outputs and saved predictions go to `rerun_outputs/`, never
to the included frozen results.

The run can take hours depending on hardware. For a numerical smoke test,
use `--replications 1`. This does not reproduce the publication tables.
`RI_THREADS` changes computational threads, not the statistical design.
Floating-point and pseudorandom differences across R builds or compilers may
prevent bitwise equality on other platforms.

### Monte Carlo transformations and seeds

For a pooled training/evaluation covariate vector v, let z=v-min(v).
The eight additional columns, in order, are `(v-mean(v))/sd(v)`,
`min(v)+max(v)-v`, `1000*v+250`, `log1p(z)`, `rank(v, ties.method="average")`,
`z^2`, `exp(v)`, and `z^3`. All original p columns remain. The two
representations append these columns to X1 and X2 respectively.

The primary experiment uses data seed `20260811+10000+r` and forest seed
`20260811+30000+r`, r=1,...,100. The seed-control experiment uses data seed
`20260812+10000+r`, baseline forest seed `20260812+30000+r`, and alternative
forest seed `20260812+50000+r`. The released code is the definitive record
of draw order. Both forest comparisons use oracle nuisances. Class sampling
certifies the pooled covariate array without outcomes, preserves ties,
canonicalizes orientation, and sorts classes by their certificates. Its
mtry is the package default applied to the number of classes. This is a
finite-array guarantee, not an out-of-sample extrapolation guarantee.

## 3. Refit the empirical comparisons

Obtain `18515_Data_and_Programs.zip` from the official
[Econometric Society page](https://www.econometricsociety.org/publications/econometrica/2022/05/01/automatic-debiased-machine-learning-causal-and-structural).
Expected archive SHA-256:
`89b0e810319eec93bc40418ac802820d27c96d3c1bf6ea9786f6b75fffa169db`.
Extract it and locate `rrr_lasso_NSW_blackbox`. The archive, its programs,
and individual-level records are not redistributed here.

Create a separate environment using
`replication/code/auto_dml/environment.yml`: R 4.0.5, randomForest 4.6-14,
nnet 7.3-16. The R package reports its version as 4.6.14.

```sh
python3 scripts/reproduce.py --stage empirical \
  --empirical-rscript /path/to/R-4.0.5/bin/Rscript \
  --package-root /path/to/rrr_lasso_NSW_blackbox
```

This rebuilds the trimmed samples, five folds, and Riesz/Lasso estimates
using the official routines. It then runs the native full and two reduced
dictionaries at seed 1, native full-dictionary seeds 1--5, and the paired
canonical full/reduced dictionaries at seed 1. All use 1,000 trees. The
canonical experiment reuses exactly the native run's cached folds and Riesz
coefficients. Numerical mtry is 8 for specification 1 and 10 for specification
2 in both canonical dictionaries. Details, rank conventions, and retained
column maps are in [EMPIRICAL.md](EMPIRICAL.md).

The complete empirical entry command was checked from newly extracted
official inputs. Rebuilt sample arrays, five folds, and Riesz coefficients
were identical to the reference caches. All 18 native fits, 30 seed-control
fits, 12 canonical fits, and their diagnostics reproduced the included
statistics exactly. Verification scope is recorded in
`metadata/verification.json`; reference cache hashes are in
`metadata/input_cache_hashes.json`. The original fold routine can warn when
the sample count is not divisible by five. Its assignment is preserved.

The reconstructed native baselines approximate the published table values,
rather than reproducing every printed decimal, because the released
original driver is configured for a different learner/specification/tuning
combination.

## 4. Build the paper

With a normal LaTeX installation containing latexmk, natbib, and booktabs:

```sh
./scripts/build_paper.sh
```

The current journal-neutral source uses 12-point type, one-inch margins,
and proofs in the appendix. The committed PDF is the delivered reference
artifact. Small line-break differences between TeX engines are possible.
