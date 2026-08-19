#!/usr/bin/env python3
"""Reconstruct the numerical claims in Tables I--III from included CSVs."""

from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def read_csv(path: str) -> list[dict[str, str]]:
    with (ROOT / path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def assert_close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if not math.isfinite(actual) or abs(actual - expected) > tolerance:
        raise AssertionError(
            f"{label}: expected {expected}, found {actual} (tol={tolerance})"
        )


def mean_se(values: list[float]) -> tuple[float, float]:
    if len(values) < 2:
        raise AssertionError("At least two replications are required")
    return statistics.mean(values), statistics.stdev(values) / math.sqrt(len(values))


def grouped_metric(
    rows: list[dict[str, str]],
    key_fields: tuple[str, ...],
) -> dict[tuple[str, ...], tuple[float, float, int]]:
    groups: dict[tuple[str, ...], list[float]] = defaultdict(list)
    for row in rows:
        key = tuple(row[field] for field in key_fields)
        groups[key].append(float(row["value"]))
    output = {}
    for key, values in groups.items():
        mean, se = mean_se(values)
        output[key] = (mean, se, len(values))
    return output


def audit_table_ii_kenya() -> None:
    rows = read_csv("results/kenya/representation_sensitivity.csv")
    by_name = {row["scenario"]: row for row in rows}
    expected = {
        "clone_hh_size_002": (17, 0.984, 5.25),
        "clone_hh_size_032": (47, 0.902, 12.76),
        "clone_hh_size_128": (143, 0.751, 21.32),
        "clone_female_head_128": (143, 0.693, 25.18),
        "clone_employed_128": (143, 0.824, 18.38),
        "clone_protein_meals_128": (143, 0.711, 23.17),
        "routine_raw_log_rank": (22, 0.978, 6.55),
        "balanced_clone_008": (128, 0.980, 6.67),
    }
    for name, (columns, rank_correlation, switching_percent) in expected.items():
        row = by_name[name]
        if int(row["p"]) != columns:
            raise AssertionError(f"Table II {name}: column count mismatch")
        assert_close(
            float(row["cate_rank_correlation"]),
            rank_correlation,
            0.00051,
            f"Table II {name} rank correlation",
        )
        assert_close(
            100 * float(row["assignment_switching"]),
            switching_percent,
            0.0051,
            f"Table II {name} assignment switching",
        )

    seed_rows = [
        row
        for row in rows
        if row["kind"] == "baseline" and row["scenario"] != "baseline"
    ]
    maximum = max(seed_rows, key=lambda row: float(row["assignment_switching"]))
    assert_close(
        100 * float(maximum["assignment_switching"]),
        0.55,
        0.0051,
        "Table II maximum seed switching",
    )
    assert_close(
        float(maximum["cate_rank_correlation"]),
        0.9999,
        0.000051,
        "Table II maximum seed rank correlation",
    )

    opportunity = read_csv("results/verification/candidate_opportunity.csv")
    by_multiplicity = {int(row["multiplicity"]): row for row in opportunity}
    assert_close(
        float(by_multiplicity[1]["expected_candidate_count"]),
        14.41,
        0.0051,
        "canonical expected candidate count",
    )
    assert_close(
        float(by_multiplicity[1]["ordinary_family_opportunity"]),
        0.9008,
        0.000051,
        "canonical singleton availability",
    )
    assert_close(
        float(by_multiplicity[128]["ordinary_family_opportunity"]),
        0.1119,
        0.000051,
        "128-copy singleton availability",
    )


def audit_tables_i_and_iii() -> None:
    robust_rows = read_csv(
        "results/monte_carlo/sign_disagreement/robustness_metrics.csv"
    )
    robust = grouped_metric(
        robust_rows,
        ("panel", "multiplicity", "trees", "method", "metric"),
    )
    metrics = (
        "sign_reversal",
        "strong_reversal_010",
        "reversal_given_true_margin_025",
    )
    panel_a = {
        1: ((17.69, 1.42), (4.00, 0.85), (1.06, 0.33)),
        2: ((24.50, 1.15), (13.38, 1.64), (3.72, 1.50)),
        4: ((30.06, 1.59), (21.00, 1.02), (8.56, 2.17)),
        8: ((35.54, 1.79), (25.54, 1.25), (15.28, 2.41)),
        16: ((41.33, 1.65), (29.50, 1.58), (22.76, 2.20)),
    }
    for multiplicity, expected_cells in panel_a.items():
        for metric, (expected_mean, expected_se) in zip(metrics, expected_cells):
            mean, se, n = robust[
                ("dose_response", str(multiplicity), "2000", "ordinary", metric)
            ]
            if n != 20:
                raise AssertionError("Table I Panel B must have 20 replications")
            assert_close(100 * mean, expected_mean, 0.0051, f"Table I Panel B m={multiplicity} {metric} mean")
            assert_close(100 * se, expected_se, 0.0051, f"Table I Panel B m={multiplicity} {metric} se")

    tree_expected = {800: 25.67, 5000: 25.76}
    for trees, expected in tree_expected.items():
        mean, _, n = robust[
            ("tree_robustness", "8", str(trees), "ordinary", "strong_reversal_010")
        ]
        if n != 20:
            raise AssertionError("Tree robustness must have 20 replications")
        assert_close(100 * mean, expected, 0.0051, f"tree robustness {trees}")

    dimension_rows = read_csv(
        "results/monte_carlo/dimension_scaling/dimension_metrics.csv"
    )
    dimension = grouped_metric(
        dimension_rows,
        ("dimension", "method", "metric"),
    )
    panel_b = {
        40: ((6.83, 1.24), (0.02, 0.01), (0.40, 0.15)),
        100: ((21.86, 1.30), (10.57, 1.82), (2.50, 0.75)),
        250: ((37.23, 1.98), (25.96, 1.51), (17.26, 2.77)),
        500: ((46.22, 0.99), (31.37, 2.10), (28.99, 1.33)),
    }
    for dimension_value, expected_cells in panel_b.items():
        for metric, (expected_mean, expected_se) in zip(metrics, expected_cells):
            mean, se, n = dimension[
                (str(dimension_value), "ordinary_default", metric)
            ]
            if n != 20:
                raise AssertionError("Table I Panel A must have 20 replications")
            assert_close(100 * mean, expected_mean, 0.0051, f"Table I Panel A p={dimension_value} {metric} mean")
            assert_close(100 * se, expected_se, 0.0051, f"Table I Panel A p={dimension_value} {metric} se")

    seed_rows = read_csv(
        "results/monte_carlo/seed_sensitivity/seed_benchmark_metrics.csv"
    )
    seed = grouped_metric(seed_rows, ("comparison", "metric"))
    panel_c = {
        "representation": ((46.10, 0.98), (31.39, 2.10), (28.83, 1.32)),
        "second_seed": ((0.09, 0.02), (0.00, 0.00), (0.02, 0.00)),
    }
    for comparison, expected_cells in panel_c.items():
        for metric, (expected_mean, expected_se) in zip(metrics, expected_cells):
            mean, se, n = seed[(comparison, metric)]
            if n != 20:
                raise AssertionError("Table I Panel C must have 20 replications")
            assert_close(100 * mean, expected_mean, 0.0051, f"Table I Panel C {comparison} {metric} mean")
            assert_close(100 * se, expected_se, 0.0051, f"Table I Panel C {comparison} {metric} se")

    for metric in metrics:
        mean, se, n = dimension[("500", "class_sampled_default", metric)]
        if n != 20 or mean != 0 or se != 0:
            raise AssertionError(f"Table I class-sampled row is not exactly zero for {metric}")


def audit_table_iii_class_sampling() -> None:
    real_rows = read_csv("results/kenya/class_sampling.csv")
    scenarios = {
        "administrative_aliases",
        "unit_recodings",
        "complements",
        "raw_log_rank",
        "routine_bundle",
    }
    ordinary_means = []
    for scenario in sorted(scenarios):
        ordinary = [
            float(row["assignment_switching"])
            for row in real_rows
            if row["scenario"] == scenario and row["method"] == "ordinary"
        ]
        quotient = [
            row
            for row in real_rows
            if row["scenario"] == scenario and row["method"] == "class_sampled"
        ]
        if len(ordinary) != 3 or len(quotient) != 3:
            raise AssertionError(f"Expected three Kenya seeds for {scenario}")
        ordinary_means.append(100 * statistics.mean(ordinary))
        for row in quotient:
            if float(row["assignment_switching"]) != 0:
                raise AssertionError(f"Class-sampled switching is nonzero for {scenario}")
            assert_close(
                float(row["cate_correlation"]),
                1.0,
                1e-14,
                f"Class-sampled CATE equality for {scenario}",
            )
    assert_close(min(ordinary_means), 5.86, 0.0051, "minimum routine Kenya switching")
    assert_close(max(ordinary_means), 7.54, 0.0051, "maximum routine Kenya switching")

    simulation_rows = read_csv(
        "results/class_sampling/simulation_natural_recodings.csv"
    )
    quotient_rows = [
        row for row in simulation_rows if row["method"] == "class_sampled"
    ]
    if not quotient_rows:
        raise AssertionError("No class-sampled simulation rows found")
    for row in quotient_rows:
        if float(row["assignment_switching"]) != 0:
            raise AssertionError("Class-sampled simulation switching is nonzero")
        if float(row["maximum_absolute_prediction_difference"]) != 0:
            raise AssertionError("Class-sampled simulation predictions differ")
        if float(row["policy_value_difference"]) != 0:
            raise AssertionError("Class-sampled simulation value difference is nonzero")

    equality_rows = read_csv("results/verification/class_sampled_reported_equality.csv")
    if len(equality_rows) != 1:
        raise AssertionError("Expected one finite-array equality row")
    equality = equality_rows[0]
    expected_design = {
        "training_observations": 400,
        "target_observations": 200,
        "trees": 400,
        "raw_dimension_canonical": 3,
        "raw_dimension_augmented": 7,
        "semantic_dimension": 3,
    }
    for field, expected in expected_design.items():
        if int(equality[field]) != expected:
            raise AssertionError(f"Table III Panel C {field} mismatch")
    for field in (
        "maximum_absolute_prediction_difference",
        "maximum_absolute_variance_difference",
        "maximum_absolute_standard_error_difference",
    ):
        if float(equality[field]) != 0:
            raise AssertionError(f"Table III Panel C {field} is nonzero")


def main() -> None:
    audit_tables_i_and_iii()
    audit_table_ii_kenya()
    audit_table_iii_class_sampling()
    print("PASS: Table I reconstructed from 20 replication-level rows per cell.")
    print("PASS: Table II reconstructed from included Kenya rows.")
    print("PASS: Table III robustness and class-sampling claims reconstructed.")
    print("PASS: candidate-opportunity calculations match the manuscript.")
    print("PASS: class-sampled predictions and assignments are exactly invariant.")
    print("PASS: all included numerical claims audited without raw Kenya data.")


if __name__ == "__main__":
    main()
