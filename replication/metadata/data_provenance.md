# Data Provenance

## Kenya cash-transfer files

The diagnostic uses household files from the public replication materials for
Egger et al. (2022), "General Equilibrium Effects of Cash Transfers:
Experimental Evidence From Kenya," DOI `10.3982/ECTA17945`. The files are
available from the original authors or journal and are mirrored at
<https://zenodo.org/records/16548593>.

The microdata are not redistributed. Exact filenames, MD5 checksums, sample
restrictions, variables, and placement instructions appear in
`replication/data/README.md`. The reproduction scripts enforce those checksums
unless the user explicitly opts into a non-replication run.

All recodings are listed in `transformations.csv`. Observations, nuisance
estimates, weights, clusters, tree counts, target share, and forest seed are
paired across each representation comparison.

The included Kenya CSVs are aggregate outputs from the checksum-verified data
run used for this repository. The raw files and fitted household-level objects
are not included. The source hashes for that run are listed in
`kenya_source_provenance.csv`; the portable scripts implement the same designs
with repository-relative paths and explicit input checksums.

## Simulation data

Simulation data are generated inside the R scripts from the schedules in
`seeds.csv`. No external inputs are required. Primary comparisons reuse the
same generated data and forest seed across paired representations and use an
independent evaluation sample. The second-seed benchmark changes only the
forest randomization.

The included Monte Carlo CSVs were regenerated on August 19, 2026, from the
public R files at their full documented settings. Run metadata and R session
information accompany each result block.

## Scope

The Kenya exercise measures representation sensitivity in a hypothetical
ranking. It does not estimate winning first resolutions, population
resolution hazards, policy accuracy, or policy welfare.
