# Simulation Provenance

The included Monte Carlo outputs were generated on August 19, 2026, from the
public R files in `replication/code/r`. The production run used R 4.3.3,
`grf` 2.4.0, the default sample sizes and tree counts declared in each script,
and twenty replications per design cell.

The run covered the dimension-scaling panel, the sign-disagreement robustness
matrix, the second-seed benchmark, and the class-sampling simulation. The
class-sampled fits had zero maximum prediction difference and zero assignment
switching across every certified recoding in the simulation.

Each result directory contains `run_metadata.txt` and `session_info.txt`.
Generated files under `rerun_outputs/` are working artifacts and are not part of
the versioned release.

`simulation_class_maps.csv` records class assignments and representatives. The
full weak-order certificate strings are omitted because they can be
reconstructed from the declared simulation array and would otherwise repeat
tens of thousands of row indices for every raw column.
