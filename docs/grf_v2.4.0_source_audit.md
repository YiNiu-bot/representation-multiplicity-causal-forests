# GRF v2.4.0 source audit

The source-faithful theorem is pinned to GRF tag `v2.4.0`, commit
`11ead50b374d3f7bd8cc88ee7205fb348efd1624`. The audit covers the following
implementation files:

- `core/src/relabeling/InstrumentalRelabelingStrategy.cpp`
- `core/src/splitting/InstrumentalSplittingRule.cpp`
- `core/src/prediction/InstrumentalPredictionStrategy.cpp`
- `core/src/prediction/collector/OptimizedPredictionCollector.cpp`
- `core/src/tree/TreeTrainer.cpp`
- `r-package/grf/bindings/CausalForestBindings.cpp`
- `r-package/grf/R/causal_forest.R`

## Translation

With unit weights, supplied `W.hat = 1/2`, supplied outcome nuisance
`mhat`, and `Z = W - 1/2`, the R wrapper passes `Y - mhat` as outcome and
`Z` as both treatment and instrument. At a node, the relabeler computes the
sample covariance ratio and sets

```text
R_i = (Z_i - mean(Z))
      * {Ytilde_i - mean(Ytilde) - theta * (Z_i - mean(Z))}.
```

The scalar `causal_forest` binding calls `instrumental_trainer`, which
instantiates `InstrumentalRelabelingStrategy`. The latter writes the response
above directly. It must not be confused with `MultiCausalRelabelingStrategy`,
whose `rho_weight` includes an inverse treatment-covariance matrix and which is
used by the separate multi-arm trainer. No factor
`1 / (n_C * Vhat_C)` enters the scalar causal-forest response.

The splitting rule maximizes `sum_L^2 / n_L + sum_R^2 / n_R` subject to:

- at least `min.node.size` observations below and above the parent instrument
  mean in each child;
- each child's instrument information being at least `alpha` times the
  parent's information;
- the optional imbalance penalty, set to zero in the theorem.

For binary randomized treatment, below and above the parent mean correspond
to controls and treated observations whenever both arms occur. Instrument
information equals `n * (1/4 - mean(Z)^2)`.

The optimized prediction collector first averages the seven leaf sufficient
statistics over nonempty trees and then applies the instrumental covariance
ratio. It does not average tree-level treatment-effect ratios. With unit
weights this is exactly equation `source-prediction` in `main.tex`.

Candidate labels are redrawn at every node. The count is Poisson with mean
`mtry`, truncated below at one and above at the number of admissible columns;
labels are then sampled without replacement.

## Exact identities

`verify_source_grf_transfer.js` checks 10,000 deterministic random arrays for:

1. the source-relabeling expansion in Appendix equation
   `app-source-response-expansion`;
2. equality of the normalized source decrease and the CART decrease of
   `4 R_i`;
3. the exact covariance-ratio prediction remainder in equation
   `source-prediction-remainder`;
4. the binary treatment-information identity.

Run with:

```sh
node replication/code/javascript/verify_source_grf_transfer.js
```

This audit verifies the software-to-notation translation. The asymptotic
probability bounds and target conclusions remain mathematical results of the
paper rather than software tests.
