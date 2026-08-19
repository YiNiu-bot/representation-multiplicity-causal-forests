# Kenya Data

The Kenya microdata are public but are not redistributed in this repository.
Obtain the replication files for Egger et al. (2022), "General Equilibrium
Effects of Cash Transfers: Experimental Evidence From Kenya," from the
authors, the journal, or the public mirror at
<https://zenodo.org/records/16548593>.

Place these two files in `replication/data/kenya/`:

| File | MD5 |
|---|---|
| `GE_HH-Analysis_AllHHs.dta` | `6841329ee31e3fee6413651bf8a1c5ff` |
| `GE_HH-Survey-BL_Analysis_AllHHs.dta` | `35b0b1384ea3cf82e38bd77c73442f3c` |

The reproduction scripts stop if a file is missing or its checksum differs.
The files are merged by `hhid_key`. The analysis keeps eligible, baselined
households with observed treatment, village, positive endline weight, and
observed endline per-capita consumption. The resulting audit sample has 4,765
households in 653 villages.

The outcome is `p2_consumption_pc_wins_PPP`, treatment is `treat`, weights are
`hhweight_EL`, villages define clusters, and `hhid_key` breaks ranking ties.
The sixteen baseline predictors and all constructed recodings are defined in
`replication/code/r/run_fixed_budget_study.R` and
`replication/code/r/run_class_sampling_study.R`.
