# Representation Multiplicity in Causal Forests

Paper and replication materials for Yi Niu, *Representation Multiplicity in
Causal Forests*. This version includes the 100-replication Monte Carlo,
CATE error and treatment-rule regret, and six Auto-DML empirical comparisons
using both native values and controlled canonical ranks.

- [Paper](paper/representation_multiplicity_causal_forests.pdf)
- [LaTeX source](paper/source/main.tex)
- [Replication instructions](replication/README.md)
- [Empirical protocol and retained columns](replication/EMPIRICAL.md)

## Verify the included results

From the repository root, with Python 3.10 or later:

```sh
python3 scripts/verify_replication.py
```

The command checks file integrity, reconstructs all four tables from
replication-level outputs, and checks paired comparisons and certificates.
It does not refit forests or certify the proofs. Full refit commands and
software requirements are in the replication instructions.

## Findings

Under the stated continuation and score-transfer conditions, the paper
characterizes representation-specific probability limits of sample-grown
honest forests. An additional failure-path depth condition gives the
intermediate attenuation formula. The CATE itself is unchanged.

Class sampling makes certified representations yield identical predictions.
In the reported Monte Carlo it reduces treatment-rule regret but increases
CATE mean squared error. The empirical aggregate ATET conclusions remain
stable, although individual first-stage contrasts change. Canonical-rank
comparisons remove nonlinear midpoint routing as an explanation, without
identifying every implementation effect or demonstrating empirical accuracy
gains.

Older versions remain in Git history. The current tree contains only the
current paper and its replication materials. Original third-party data and
replication programs must be obtained from their distributors.

See [rights](RIGHTS.md), [third-party notices](THIRD_PARTY_NOTICES.md), and
[citation metadata](CITATION.cff).
