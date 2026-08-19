# Representation Multiplicity in Causal Forests

This repository contains the paper, source files, numerical outputs,
and replication code for **Representation Multiplicity in Causal Forests**.

- [Paper PDF](paper/representation_multiplicity_causal_forests.pdf)
- [LaTeX source](paper/source/main.tex)
- [Replication guide](replication/README.md)
- [`grf` 2.4.0 source audit](docs/grf_v2.4.0_source_audit.md)

## Main result

Feature-subsampled forests sample covariate columns rather than economic
information. Retaining several split-equivalent encodings can therefore change
which split actions are available. The paper characterizes the resulting
causal-forest targets, gives a signed example with representation-dependent
limits, and studies a class-sampled correction on a declared finite array.

The numerical exercises establish sensitivity in the stated designs. They do
not estimate how often the problem occurs in applied work, and the Kenya
ranking exercise is not a welfare analysis.

## Verify the repository

The fast verifier requires Python 3 and Node.js. It checks every versioned file,
reconstructs Tables I--III from replication-level CSVs, and runs the
deterministic proof and implementation checks.

```sh
./scripts/verify_release.sh
```

The R smoke tests additionally require R 4.3.3 and `grf` 2.4.0:

```sh
export RI_RLIB="$PWD/.r-library"
Rscript scripts/install_r_dependencies.R
./scripts/smoke_test.sh
```

The included Monte Carlo outputs were regenerated with the public code. Full
Monte Carlo and Kenya reruns remain separate because they are substantially
more expensive and the Kenya microdata are not redistributed. See the
[replication guide](replication/README.md) for exact commands.

## Repository layout

```text
paper/          PDF and journal-neutral LaTeX source
replication/    Code, included outputs, data instructions, and metadata
scripts/        Build, verification, and reproduction entry points
docs/           Source audit for the implementation-specific theorem
```

The raw Kenya microdata are not redistributed. Acquisition instructions and
input checksums are provided in [`replication/data/README.md`](replication/data/README.md).

## Citation and rights

Citation metadata are in [`CITATION.cff`](CITATION.cff). Author-written code is
released under the MIT License in [`LICENSE-CODE`](LICENSE-CODE). The paper,
documentation, included results, and third-party inputs have separate rights
described in [`RIGHTS.md`](RIGHTS.md) and
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
