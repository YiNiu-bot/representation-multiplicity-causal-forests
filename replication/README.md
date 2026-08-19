# Replication Guide

## What is included

`results/` contains the replication-level CSVs used in Tables I--III. These
files are sufficient to reconstruct every displayed number without the Kenya
microdata. `metadata/claims_to_outputs.csv` maps each table or numerical claim
to its included output, generating code, and audit.

The Monte Carlo outputs were regenerated on August 19, 2026, with the public R
files in this repository, R 4.3.3, and `grf` 2.4.0. Their run metadata and R
session information are stored beside the CSVs. The Kenya aggregate outputs
come from a checksum-verified restricted-data run because the microdata are not
redistributed. Hashes of the source files used for that run are recorded in
`metadata/kenya_source_provenance.csv`.

The `grf` source audit is pinned to tag `v2.4.0`, commit
`11ead50b374d3f7bd8cc88ee7205fb348efd1624`.

The public R files use repository-relative paths and documented environment
variables. A different R build, compiler, or numerical library can produce
small stochastic-fit differences even with the same package version and seeds.
The repository verifier therefore checks the included replication-level outputs,
while a full refit is a separate computational audit.

## Verification levels

### 1. Included-output integrity

From the repository root:

```sh
./scripts/verify_release.sh
```

This fast check verifies file hashes, reconstructs Tables I--III, and runs the
deterministic JavaScript checks. It does not refit forests.

### 2. R implementation smoke tests

```sh
export RI_RLIB="$PWD/.r-library"
Rscript scripts/install_r_dependencies.R
./scripts/smoke_test.sh
```

This parses every R file, tests the split-class implementation, and reruns the
finite-array prediction and standard-error equality check.

### 3. Self-contained Monte Carlo reruns

```sh
./scripts/reproduce_simulations.sh
```

The full command is computationally expensive. It writes new files under
`rerun_outputs/` and never overwrites the included results. Sample sizes,
replication counts, tree counts, and threads can be reduced with the `RI_*`
environment variables defined near the top of each R script.

### 4. Kenya reruns

Follow `data/README.md`, then run:

```sh
./scripts/reproduce_kenya.sh replication/data/kenya
```

This executes the fixed-budget audit and the class-sampling study. The scripts
verify both input MD5 checksums before fitting any model.

## Designs and interpretation

The Monte Carlo uses randomized treatment, paired samples and seeds, oracle
nuisance functions, and independent evaluation samples. Its primary outcomes
are paired CATE sign disagreement and margin-restricted sign disagreement.

The Kenya exercise ranks households by out-of-bag causal-forest scores and
records the fraction whose top-half selection status changes. It is a
representation-sensitivity diagnostic. Because the original experiment has
clustered saturation and interference, the stored comparisons are not
interpreted as estimates of policy welfare.

The class-sampled procedure certifies weak-order split classes jointly on the
training and declared evaluation arrays. It does not merge merely correlated
predictors and does not claim invariance at undeclared future points.
