# Empirical protocol

The six cells are NSW, PSID, and CPS comparison samples, each under
specifications 1 and 2 of Chernozhukov, Newey, and Singh (2022),
*Automatic Debiased Machine Learning of Causal and Structural Effects*.
Only the random-forest outcome-regression dictionary is changed.

## Sample and score

Official `get_data_intersection` constructs and trims each sample. Its
outcome, treatment, trimming, and underlying covariates are held fixed.
Five folds are generated once with seed 1 by the original recorded draw
order. The Riesz/Lasso step uses the original dictionary and theoretical
penalty convention `c_rr=0.5`. Fold-specific Riesz coefficients are cached
and reused without changes for every representation and seed comparison.
The original orthogonal ATET score and standard-error formula are retained.

The outcome dictionary is `[1, W, X, W*X]`. Specification 1 uses age,
education, married, black, Hispanic, earnings in 1974 and 1975, then squared
age, education, and the two earnings variables. Specification 2 adds
indicators for zero 1974 earnings, zero 1975 earnings, and no degree.
Full dictionaries therefore have 24 and 30 columns. Constants remain in
both dictionaries and are not combined with nonconstant classes.

## Certification and retained columns

Certification pools three arrays: observed treatment and both treatment
counterfactuals for every trimmed covariate row. It uses covariates and
treatment settings, never outcomes. Equivalent columns must induce the
same weak ordering, including ties, up to reversal on this whole array.
The smallest original column index is retained in each class. Columns
remaining after reduction are 19/25 in NSW and CPS and 20/26 in PSID.

Every original member, representative, orientation, and class size is
listed in `results/empirical/native/equivalence_classes.csv` and the six
`results/empirical/canonical/*_certificate.csv` files. These files fully
specify the retained dictionaries. Observed-row-only certificates can
misclassify treatment interactions and must not replace the pooled check.

## Native comparison

Native forests use original numerical values and original midpoint routing.
At seed 1, the reduced dictionary is fitted with both the recomputed
regression default mtry (6 or 8) and the fixed original numerical value
(8 or 10). The full-dictionary seed benchmark compares seed 1 with seeds
2--5. These results remain separate from the canonical experiment.

## Canonical-rank comparison

For each class, the representative's distinct pooled values are sorted and
assigned consecutive integer ranks starting at 1. Ties receive the same
rank. Each member column, including an order-reversing member, is mapped
to this representative orientation. Thus all equivalent columns have
identical numeric values on every observed and counterfactual row.

The full dictionary retains all these columns, including duplicates. The
reduced dictionary keeps the original representative indices. Both use
identical rank maps, seed 1, five folds, 1,000 trees, and numerical mtry 8
or 10. The mapping is applied only to forest fitting and prediction.
The original raw dictionary continues to enter the Riesz score. No Riesz
coefficients, outcomes, folds, or score terms are re-estimated for this
comparison. Exact rank equality is checked for all class members.

This removes nonlinear midpoint routing between equivalent columns as an
alternative explanation. With a fixed numerical seed, changing the number
of columns can still change tie handling and consumption of pseudorandom
draws. The experiment is not a decomposition of all implementation effects.
It also does not measure prediction error or welfare, because the true
conditional effects are unknown.

## Diagnostics

The individual contrast is the cross-fitted outcome regression
`gamma(1,x)-gamma(0,x)`. It is distinct from the orthogonal ATET score.
Sign disagreement is the fraction with unequal `sign()` values,
reported over all observations and over treated observations separately.
Aggregate changes are reduced minus full ATET in dollars. Every reported
empirical fit uses R 4.0.5 and randomForest 4.6-14. Detailed protocols and
session information accompany the included results.
