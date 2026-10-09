#!/usr/bin/env python3
"""Reconstruct both manuscript tables from the released CSVs."""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "replication/metadata/manifest.json"


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def numeric_cells(line):
    content = line.split(r"\\", 1)[0]
    return [float(re.sub(r"[^0-9.+-]", "", cell))
            for cell in content.split(" & ") if re.search(r"\d", cell)]


def mean_and_mcse(values):
    require(len(values) == 10, "Every Monte Carlo cell must have ten replications.")
    return statistics.mean(values), statistics.stdev(values) / math.sqrt(len(values))


def check_table_one(source):
    metrics = read_csv(ROOT / "replication/results/monte_carlo/metrics.csv")
    comparisons = read_csv(ROOT / "replication/results/monte_carlo/comparisons.csv")
    grouped = defaultdict(dict)
    paired = defaultdict(dict)
    for row in metrics:
        key = (int(row["node_size"]), int(row["multiplicity"]), row["variant"], row["metric"])
        replication = int(row["replication"])
        require(replication not in grouped[key], "Duplicate Monte Carlo observation.")
        grouped[key][replication] = float(row["value"])
    for row in comparisons:
        if row["metric"] == "sign_disagreement":
            key = (int(row["node_size"]), int(row["multiplicity"]), row["regime"])
            replication = int(row["replication"])
            require(replication not in paired[key], "Duplicate paired comparison.")
            paired[key][replication] = float(row["value"])
    for values in list(grouped.values()) + list(paired.values()):
        require(set(values) == set(range(1, 11)), "Incomplete Monte Carlo replications.")

    table = source.split("% BEGIN GENERATED NODE GRID TABLE", 1)[1].split(
        "% END GENERATED NODE GRID TABLE", 1)[0].splitlines()
    rows = [i for i, line in enumerate(table) if re.match(r"^\d+ &", line)]
    require(len(rows) == 15, "Table 1 must contain fifteen rows.")
    index = 0
    for node in (5, 10, 20):
        for multiplicity in (1, 5, 10, 40, 400):
            means, errors = [], []
            for regime in ("fixed", "default", "linear"):
                cells = [paired[node, multiplicity, regime]]
                cells.extend(grouped[node, multiplicity, f"{regime}_x{h}", metric]
                             for metric in ("mse", "regret") for h in (1, 2))
                for j, values in enumerate(cells):
                    mean, mcse = mean_and_mcse([100 * value for value in values.values()])
                    digits = 1 if j == 0 else 2
                    means.append(round(mean, digits))
                    errors.append(round(mcse, digits))
                for metric in ("mse", "regret", "wrong_sign"):
                    for replication in range(1, 11):
                        reference = grouped[node, 1, f"class_{regime}", metric][replication]
                        current = grouped[node, multiplicity, f"class_{regime}", metric][replication]
                        require(current == reference, "Class-sampled results must coincide.")
                if multiplicity == 1:
                    require(all(value == 0 for value in cells[0].values()), "Baseline SR must be zero.")
            row = rows[index]
            require(numeric_cells(table[row]) == [float(multiplicity), *means],
                    f"Table 1 means differ at node size {node}, multiplicity {multiplicity}.")
            require(numeric_cells(table[row + 1]) == errors,
                    f"Table 1 MCSEs differ at node size {node}, multiplicity {multiplicity}.")
            index += 1


def check_table_two(source):
    results = read_csv(ROOT / "replication/results/empirical/results.csv")
    seeds = read_csv(ROOT / "replication/results/empirical/seed_results.csv")
    require(len(results) == 9 and len(seeds) == 36, "Expected nine applications and four seeds each.")
    table = source.split("% BEGIN GENERATED EMPIRICAL TABLE", 1)[1].split(
        "% END GENERATED EMPIRICAL TABLE", 1)[0].splitlines()
    require(r"\(|\Delta\mathrm{ATET}|\)" in "\n".join(table), "Seed column must show absolute ATET changes.")
    start = next(i for i, line in enumerate(table) if line.startswith("NSW & 1"))
    lines = [line for line in table[start:] if " & " in line]
    require(len(lines) == 27, "Table 2 must contain three lines for each application.")
    for index, row in enumerate(results):
        dataset, spec = row["dataset"], int(row["specification"])
        selected = [seed for seed in seeds if seed["dataset"] == dataset
                    and int(seed["specification"]) == spec]
        require({int(seed["seed"]) for seed in selected} == {2, 3, 4, 5}, "Incomplete seed benchmark.")
        full = float(row["full_atet"])
        deltas = [float(seed["abs_atet_change"]) for seed in selected]
        sign_changes = [float(seed["sr"]) for seed in selected]
        for field, expected in (("seed_delta_mean", statistics.mean(deltas)),
                                ("seed_delta_max", max(deltas)),
                                ("seed_sr_mean", statistics.mean(sign_changes)),
                                ("seed_sr_max", max(sign_changes))):
            require(abs(float(row[field]) - expected) < 1e-8, f"Seed summary differs for {dataset} {spec}.")
        expected_first = [float(spec), round(full)]
        expected_second = [round(float(row["full_se"]))]
        expected_third = []
        for method in ("default", "fixed", "class"):
            atet, delta = float(row[f"{method}_atet"]), float(row[f"{method}_delta"])
            require(abs(delta - (atet - full)) < 1e-8, "ATET change differs from the level comparison.")
            expected_first.extend((round(atet), round(float(row[f"{method}_sr"]), 1)))
            expected_second.append(round(float(row[f"{method}_se"])))
            expected_third.append(round(delta))
            significant = abs(atet / float(row[f"{method}_se"])) > 1.959963984540054
            baseline = abs(full / float(row["full_se"])) > 1.959963984540054
            require(atet > 0 and significant == baseline, "ATET significance classification changed.")
        expected_first.extend((round(float(row["seed_delta_mean"])), round(float(row["seed_sr_mean"]), 1)))
        expected_second.extend((round(float(row["seed_delta_max"])), round(float(row["seed_sr_max"]), 1)))
        for offset, expected in enumerate((expected_first, expected_second, expected_third)):
            require(numeric_cells(lines[3 * index + offset]) == expected,
                    f"Table 2 differs for {dataset} specification {spec}, line {offset + 1}.")
        require(float(row["class_sr_between_representations"]) == 0, "Class-sampled SR must be zero.")
        require(f"{float(row['full_se']):.2f}" == f"{float(row['published_se']):.2f}", "Published SE mismatch.")
        if (dataset, spec) == ("CPS", 3):
            require(f"{full:.2f}" == "1906.63", "CPS specification 3 baseline changed.")
        else:
            require(f"{full:.2f}" == f"{float(row['published_atet']):.2f}", "Published ATET mismatch.")


def payload_hashes():
    paths = sorted(path for folder in (ROOT / "paper", ROOT / "replication", ROOT / "scripts")
                   for path in folder.rglob("*")
                   if path.is_file() and path != MANIFEST
                   and path != ROOT / "paper/source/main.pdf"
                   and path.suffix in {".csv", ".R", ".py", ".sh", ".tex", ".pdf", ".md"})
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    source = (ROOT / "paper/source/main.tex").read_text()
    check_table_one(source)
    check_table_two(source)
    for row in read_csv(ROOT / "replication/metadata/empirical_checks.csv"):
        require(row["passed"] == "TRUE", "A saved empirical fit check failed.")
    hashes = payload_hashes()
    if args.write_manifest:
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        MANIFEST.write_text(json.dumps({"date": "2026-10-09", "sha256": hashes}, indent=2) + "\n")
    else:
        manifest = json.loads(MANIFEST.read_text())
        require(manifest["sha256"] == hashes, "Replication file hashes differ from the manifest.")
    print("Verified Table 1 (15 rows), Table 2 (9 rows), class invariance, seed benchmarks, and file hashes.")


if __name__ == "__main__":
    main()
