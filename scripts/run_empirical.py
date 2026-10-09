#!/usr/bin/env python3
"""Run the original and paired empirical fits."""

import argparse
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("rerun_outputs/empirical"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    code = root / "replication/code/auto_dml"
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    paired = output / "paired_coordinates"
    spec3 = output / "spec3_known_multiplicity"

    def run(name, *arguments):
        subprocess.run([args.rscript, str(code / name), *map(str, arguments)], check=True)

    run("run_literal_rf.R", args.package_root.resolve(), output)
    run("run_paired_coordinates.R", output, paired)
    run("verify_paired_coordinates.R", paired)
    run("run_spec3_known_multiplicity.R", output, spec3)
    run("verify_spec3_known_multiplicity.R", output, spec3)


if __name__ == "__main__":
    main()
