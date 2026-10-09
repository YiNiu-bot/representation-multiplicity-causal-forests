#!/usr/bin/env python3
"""Run the multiplicity, mtry, and node-size grid."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--rlib")
    parser.add_argument("--output", type=Path, default=Path("rerun_outputs/monte_carlo"))
    parser.add_argument("--replications", type=int, default=10)
    parser.add_argument("--multiplicities", default="5,10,40,400")
    parser.add_argument("--node-sizes", default="5,10,20")
    parser.add_argument("--trees", type=int, default=2000)
    parser.add_argument("--n-train", type=int, default=3000)
    parser.add_argument("--n-test", type=int, default=10000)
    parser.add_argument("--workers", type=int, default=max((os.cpu_count() or 1) - 1, 1))
    args = parser.parse_args()
    if args.replications < 3 or args.workers < 1:
        parser.error("Use at least three replications and one worker.")
    nodes = [int(value) for value in args.node_sizes.split(",")]
    multiplicities = [int(value) for value in args.multiplicities.split(",")]
    if min(nodes + multiplicities) < 1:
        parser.error("Node sizes and multiplicities must be positive.")
    root = Path(__file__).resolve().parents[1]
    script = root / "replication/code/monte_carlo/run_multiplicity_node_grid.R"
    output = args.output.resolve()
    for folder in (output, output / "references", output / "cells", output / "logs"):
        folder.mkdir(parents=True, exist_ok=True)
    common = [args.rscript, str(script), f"--output={output}",
              f"--trees={args.trees}", f"--n-train={args.n_train}", f"--n-test={args.n_test}"]
    if args.rlib:
        common.append(f"--rlib={args.rlib}")
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       VECLIB_MAXIMUM_THREADS="1")

    def run(job):
        name, options = job
        with (output / "logs" / f"{name}.log").open("w") as log:
            subprocess.run(common + options, stdout=log, stderr=subprocess.STDOUT,
                           env=environment, check=True)
        print(f"Completed {name}", flush=True)

    run(("prepare", ["--mode=prepare"]))
    references = [(f"reference_k{k}_r{r}", ["--mode=reference", f"--node-size={k}",
                   f"--replication={r}"])
                  for k in nodes for r in range(1, args.replications + 1)]
    cells = [(f"cell_k{k}_c{c}_r{r}", ["--mode=worker", f"--node-size={k}",
              f"--multiplicity={c}", f"--replication={r}", "--verify-classes=1"])
             for k in nodes for c in multiplicities for r in range(1, args.replications + 1)]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(run, references))
        list(pool.map(run, cells))
    run(("summary", ["--mode=summarize", f"--replications={args.replications}",
                     f"--multiplicities={args.multiplicities}", f"--node-sizes={args.node_sizes}"]))
    print(f"Results saved in {output}")


if __name__ == "__main__":
    main()
