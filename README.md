# Representation Multiplicity in Causal Forests

Paper and replication materials for Yi Niu, *Representation Multiplicity in
Causal Forests*.

- [Paper](paper/representation_multiplicity_causal_forests.pdf)
- [LaTeX source](paper/source/main.tex)
- [Replication instructions](replication/README.md)

The Monte Carlo varies representation multiplicity, `mtry`, and minimum node
size while holding six underlying covariates fixed. The application replicates
specifications 1--3 of Chernozhukov, Newey, and Singh (2022) using NSW, PSID,
and CPS. The repository includes the results behind both tables.

Run the table and file checks with Python 3.10 or later.

```sh
python3 scripts/verify_replication.py
```

Author-written code is covered by [LICENSE-CODE](LICENSE-CODE). Original
empirical data and programs are obtained separately as described in the
replication instructions.
