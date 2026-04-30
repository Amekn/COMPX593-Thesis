#!/usr/bin/env python3
"""Validate thesis concordance tables and seeded figure inputs.

This script is intentionally focused on the highest-risk metric path used by
figures 31-34 and F3-F6. It checks that the notebook helper, the compiled
VariantConcordance binary, cached CSV tables, and seeded figure inputs agree on
the threshold semantics and basic numerical invariants.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


THRESHOLD_METRIC_COLUMNS = (
    "source_variants",
    "groundtruth_variants",
    "intersection_variants",
    "union_variants",
    "exact_overlap_mass",
    "weighted_jaccard",
    "jensen_shannon_similarity",
)

DERIVED_THRESHOLD_COLUMNS = ("precision", "recall", "f1")

SEEDED_COMPARE_COLUMNS = (
    "weighted_jaccard",
    "jensen_shannon_similarity",
    "exact_overlap_mass",
    "intersection_variants",
    "source_variants",
    "groundtruth_variants",
    "precision",
    "recall",
    "f1",
)

NORMALISED_AUC_COLUMNS = (
    "norm_weighted_jaccard_auc",
    "norm_jensen_shannon_similarity_auc",
    "norm_exact_overlap_mass_auc",
    "norm_f1_auc",
)

SEEDED_RUN_MANIFEST_COLUMNS = (
    "name",
    "base_model",
    "family",
    "variant",
    "seed_index",
    "seed_label",
    "seed",
    "training_input",
    "min_accuracy_save_ctc",
    "lr",
    "batch",
    "training_chunks",
    "validation_chunks",
    "run_dir",
    "model_wrapper_dir",
)

SEEDED_SEED_MANIFEST_COLUMNS = ("seed_index", "seed_label", "seed")


class ValidationFailure(RuntimeError):
    """Raised when a validation check fails."""


def find_repo_root(start: Path) -> Path:
    for candidate in (start.resolve(), *start.resolve().parents):
        if (candidate / "tables" / "model.csv").exists() and (
            candidate / "src" / "ipynb" / "plot.ipynb"
        ).exists():
            return candidate
    raise ValidationFailure(f"Could not locate repository root from {start}")


def ok(message: str) -> None:
    print(f"[ok] {message}")


def skip(message: str) -> None:
    print(f"[skip] {message}")


def fail(message: str) -> None:
    raise ValidationFailure(message)


def read_notebook(repo_root: Path) -> dict[str, Any]:
    notebook_path = repo_root / "src" / "ipynb" / "plot.ipynb"
    return json.loads(notebook_path.read_text())


def notebook_code_cells(notebook: dict[str, Any]) -> list[tuple[int, str]]:
    cells = []
    for index, cell in enumerate(notebook["cells"]):
        if cell.get("cell_type") == "code":
            cells.append((index, "".join(cell.get("source", []))))
    return cells


def validate_notebook_syntax(repo_root: Path) -> None:
    notebook = read_notebook(repo_root)
    for index, source in notebook_code_cells(notebook):
        ast.parse(source, filename=f"plot.ipynb:cell-{index}")
    ok("plot.ipynb code cells parse as Python")


def load_plot_namespace(repo_root: Path) -> dict[str, Any]:
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.chdir(repo_root)
    notebook = read_notebook(repo_root)
    namespace: dict[str, Any] = {"__name__": "__validate_plot_notebook__"}

    # Cell 2 contains imports/constants. Cell 4 contains the shared helpers and
    # top-level table dataframes used by later figure cells.
    for index in (2, 4):
        source = "".join(notebook["cells"][index]["source"])
        source = source.replace(
            "if PANEL_MANIFEST_PATH.exists():\n    PANEL_MANIFEST_PATH.unlink()",
            "if PANEL_MANIFEST_PATH.exists():\n    pass",
        )
        exec(compile(source, f"plot.ipynb:cell-{index}", "exec"), namespace)

    return namespace


def read_variant_counts(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    if "variant_key" not in df or "count" not in df:
        fail(f"{path} is missing variant_key/count columns")
    return df.set_index("variant_key")["count"].astype(float).sort_values(ascending=False)


def compare_float(actual: float, expected: float, *, label: str, tolerance: float = 1e-9) -> None:
    if not (math.isfinite(actual) and math.isfinite(expected)):
        if math.isnan(actual) and math.isnan(expected):
            return
        fail(f"{label}: non-finite mismatch actual={actual!r}, expected={expected!r}")
    if abs(actual - expected) > tolerance:
        fail(f"{label}: actual={actual:.12g}, expected={expected:.12g}, diff={abs(actual - expected):.3g}")


def expected_threshold_semantics_metrics() -> dict[str, float]:
    source = np.array([10.0, 2.0, 0.0])
    truth = np.array([5.0, 1.0, 2.0])
    source_prob = source / source.sum()
    truth_prob = truth / truth.sum()
    midpoint = 0.5 * (source_prob + truth_prob)

    def kl(prob: np.ndarray, ref: np.ndarray) -> float:
        mask = prob > 0
        return float(np.sum(prob[mask] * np.log2(prob[mask] / ref[mask])))

    js_divergence = 0.5 * kl(source_prob, midpoint) + 0.5 * kl(truth_prob, midpoint)
    return {
        "source_variants": 2.0,
        "groundtruth_variants": 2.0,
        "intersection_variants": 1.0,
        "union_variants": 3.0,
        "exact_overlap_mass": float(np.minimum(source_prob, truth_prob).sum()),
        "weighted_jaccard": float(np.minimum(source, truth).sum() / np.maximum(source, truth).sum()),
        "jensen_shannon_similarity": 1.0 - js_divergence,
    }


def validate_threshold_semantics(namespace: dict[str, Any], repo_root: Path) -> None:
    expected = expected_threshold_semantics_metrics()
    source_counts = pd.Series({"A": 10.0, "B": 2.0, "C": 1.0})
    truth_counts = pd.Series({"A": 5.0, "B": 1.0, "D": 2.0})
    metrics = namespace["compute_variant_concordance_from_counts"](
        source_counts,
        truth_counts,
        threshold=2,
    )
    for column, expected_value in expected.items():
        compare_float(float(metrics[column]), expected_value, label=f"notebook threshold semantics {column}")

    concordance_bin = repo_root / "build" / "VariantConcordance"
    if not concordance_bin.exists():
        fail(f"Compiled VariantConcordance binary is missing: {concordance_bin}")

    with tempfile.TemporaryDirectory(prefix="variant-concordance-semantics.") as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        source_path = tmp_dir / "source.variants.tsv"
        truth_path = tmp_dir / "truth.variants.tsv"
        source_path.write_text("variant_key\tcount\nA\t10\nB\t2\nC\t1\n")
        truth_path.write_text("variant_key\tcount\nA\t5\nB\t1\nD\t2\n")
        result = subprocess.run(
            [str(concordance_bin), str(source_path), str(truth_path), "2"],
            check=True,
            capture_output=True,
            text=True,
        )

    rows = list(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
    if len(rows) != 1:
        fail(f"VariantConcordance synthetic check returned {len(rows)} rows")
    output = rows[0]
    cpp_mapping = {
        "source_variants": "source_haplotypes",
        "groundtruth_variants": "groundtruth_haplotypes",
        "intersection_variants": "intersection_haplotypes",
        "union_variants": "union_haplotypes",
        "exact_overlap_mass": "exact_overlap_mass",
        "weighted_jaccard": "weighted_jaccard",
        "jensen_shannon_similarity": "jensen_shannon_similarity",
    }
    for metric_name, output_name in cpp_mapping.items():
        compare_float(
            float(output[output_name]),
            expected[metric_name],
            label=f"C++ threshold semantics {metric_name}",
        )
    ok("threshold opposite-side mass semantics match expected values in notebook and C++")


def validate_notebook_concordance_matches_correlation(namespace: dict[str, Any], tolerance: float) -> None:
    correlation_df = namespace["correlation_df"]
    resolve_variant_path = namespace["resolve_variant_path"]
    compute_metrics = namespace["compute_variant_concordance_from_counts"]
    truth_counts = read_variant_counts(resolve_variant_path("UDP0057"))
    max_diff = 0.0
    worst_label = ""

    for _, row in correlation_df.iterrows():
        display_name = row["display_name"]
        source_counts = read_variant_counts(resolve_variant_path(display_name))
        for threshold in range(1, 11):
            metrics = compute_metrics(source_counts, truth_counts, threshold=threshold)
            source_variants = float(metrics["source_variants"])
            groundtruth_variants = float(metrics["groundtruth_variants"])
            intersection = float(metrics["intersection_variants"])
            precision = intersection / source_variants if source_variants else math.nan
            recall = intersection / groundtruth_variants if groundtruth_variants else math.nan
            f1 = 2 * precision * recall / (precision + recall) if precision + recall else math.nan
            derived = {"precision": precision, "recall": recall, "f1": f1}

            for column in THRESHOLD_METRIC_COLUMNS:
                table_column = f"{column}_{threshold}"
                if table_column not in row or pd.isna(row[table_column]):
                    continue
                actual = float(metrics[column])
                expected = float(row[table_column])
                diff = abs(actual - expected)
                if diff > max_diff:
                    max_diff = diff
                    worst_label = f"{display_name} threshold={threshold} {column}"
            for column in DERIVED_THRESHOLD_COLUMNS:
                table_column = f"{column}_{threshold}"
                if table_column not in row or pd.isna(row[table_column]):
                    continue
                actual = float(derived[column])
                expected = float(row[table_column])
                diff = abs(actual - expected)
                if diff > max_diff:
                    max_diff = diff
                    worst_label = f"{display_name} threshold={threshold} {column}"

    if max_diff > tolerance:
        fail(f"notebook concordance helper differs from correlation.csv at {worst_label}: {max_diff:.3g}")
    ok(f"notebook concordance helper matches correlation.csv across thresholds; max abs diff {max_diff:.3g}")


def coerce_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for column in df.columns:
        if df[column].dtype == object:
            converted = pd.to_numeric(df[column], errors="coerce")
            non_null = df[column].notna()
            if bool(converted[non_null].notna().all()):
                df[column] = converted
    return df


def require_columns(df: pd.DataFrame, path: Path, columns: tuple[str, ...]) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        fail(f"{path} is missing required columns: {missing}")


def max_dataframe_difference(left: pd.DataFrame, right: pd.DataFrame, *, label: str) -> float:
    if list(left.columns) != list(right.columns):
        missing_left = sorted(set(right.columns) - set(left.columns))
        missing_right = sorted(set(left.columns) - set(right.columns))
        fail(f"{label}: column mismatch missing_left={missing_left} missing_right={missing_right}")
    if len(left) != len(right):
        fail(f"{label}: row count mismatch {len(left)} != {len(right)}")

    max_diff = 0.0
    for column in left.columns:
        left_numeric = pd.to_numeric(left[column], errors="coerce").astype(float)
        right_numeric = pd.to_numeric(right[column], errors="coerce").astype(float)
        numeric_mask = left_numeric.notna() | right_numeric.notna()
        if numeric_mask.any():
            diff = (left_numeric[numeric_mask] - right_numeric[numeric_mask]).abs()
            if diff.notna().any():
                max_diff = max(max_diff, float(diff.max()))
            both_nan = left_numeric[numeric_mask].isna() & right_numeric[numeric_mask].isna()
            if not bool(both_nan.all() or (diff.fillna(0.0) <= 1e-9).all()):
                fail(f"{label}: numeric mismatch in column {column}")
        text_mask = ~numeric_mask
        if text_mask.any():
            left_text = left.loc[text_mask, column].fillna("").astype(str)
            right_text = right.loc[text_mask, column].fillna("").astype(str)
            if not bool((left_text == right_text).all()):
                fail(f"{label}: text mismatch in column {column}")
    return max_diff


def validate_generate_correlation(repo_root: Path, tolerance: float) -> None:
    script_path = repo_root / "src" / "bash" / "GenerateCorrelationCsv.sh"
    input_path = repo_root / "tables" / "correlation.csv"
    test_root = repo_root / "src" / "ipynb" / "plot" / "test"
    concordance_bin = repo_root / "build" / "VariantConcordance"
    if not concordance_bin.exists():
        fail(f"Compiled VariantConcordance binary is missing: {concordance_bin}")

    with tempfile.NamedTemporaryFile(prefix="correlation-regenerated.", suffix=".csv", delete=False) as handle:
        output_path = Path(handle.name)
    try:
        result = subprocess.run(
            [
                str(script_path),
                str(input_path),
                str(test_root),
                str(output_path),
                "--concordance-bin",
                str(concordance_bin),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            stderr_tail = "\n".join(result.stderr.splitlines()[-20:])
            fail(f"GenerateCorrelationCsv.sh failed with exit {result.returncode}:\n{stderr_tail}")
        expected = coerce_csv(input_path)
        actual = coerce_csv(output_path)
        max_diff = max_dataframe_difference(actual, expected, label="regenerated correlation.csv")
        if max_diff > tolerance:
            fail(f"regenerated correlation.csv differs from checked-in table by {max_diff:.3g}")
        ok(f"GenerateCorrelationCsv.sh reproduces tables/correlation.csv; max abs diff {max_diff:.3g}")
    finally:
        output_path.unlink(missing_ok=True)


def validate_seeded_manifests(repo_root: Path, cached: pd.DataFrame) -> pd.DataFrame:
    seeded_tables_dir = repo_root / "tables" / "seeded"
    run_manifest_path = seeded_tables_dir / "run_manifest.csv"
    seed_manifest_path = seeded_tables_dir / "seed_manifest.csv"
    if not run_manifest_path.exists():
        fail(f"seeded run manifest is missing: {run_manifest_path}")
    if not seed_manifest_path.exists():
        fail(f"seed manifest is missing: {seed_manifest_path}")

    run_manifest = coerce_csv(run_manifest_path)
    seed_manifest = coerce_csv(seed_manifest_path)
    require_columns(run_manifest, run_manifest_path, SEEDED_RUN_MANIFEST_COLUMNS)
    require_columns(seed_manifest, seed_manifest_path, SEEDED_SEED_MANIFEST_COLUMNS)
    require_columns(cached, seeded_tables_dir / "threshold_long.csv", ("name", "seed_index", "seed_label", "seed", "threshold"))

    duplicated_runs = run_manifest.loc[run_manifest["name"].duplicated(), "name"].tolist()
    if duplicated_runs:
        fail(f"seeded run manifest contains duplicate run names: {duplicated_runs[:5]}")

    duplicated_seed_rows = seed_manifest.loc[seed_manifest["seed_index"].duplicated(), "seed_index"].tolist()
    if duplicated_seed_rows:
        fail(f"seed manifest contains duplicate seed_index values: {duplicated_seed_rows[:5]}")

    manifest_seed_rows = (
        run_manifest[["seed_index", "seed_label", "seed"]]
        .drop_duplicates()
        .sort_values("seed_index")
        .reset_index(drop=True)
    )
    seed_manifest_rows = seed_manifest.sort_values("seed_index").reset_index(drop=True)
    max_dataframe_difference(seed_manifest_rows, manifest_seed_rows, label="seed manifest vs run manifest seed rows")

    run_names = set(run_manifest["name"].astype(str))
    cached_names = set(cached["name"].astype(str))
    if cached_names != run_names:
        missing_cached = sorted(run_names - cached_names)
        missing_manifest = sorted(cached_names - run_names)
        fail(
            "seeded threshold cache/run manifest name mismatch "
            f"missing_cached={missing_cached[:5]} missing_manifest={missing_manifest[:5]}"
        )

    cached_keys = cached[["name", "threshold"]].copy()
    cached_keys["threshold"] = cached_keys["threshold"].astype(int)
    duplicated_cached = cached_keys.loc[cached_keys.duplicated(), "name"].tolist()
    if duplicated_cached:
        fail(f"seeded threshold cache contains duplicate name/threshold rows: {duplicated_cached[:5]}")

    expected_thresholds = set(range(1, 11))
    observed_thresholds = set(cached_keys["threshold"])
    if observed_thresholds != expected_thresholds:
        fail(f"seeded threshold cache thresholds are {sorted(observed_thresholds)}; expected 1..10")

    threshold_counts = cached_keys.groupby("name")["threshold"].nunique()
    if not bool((threshold_counts == len(expected_thresholds)).all()):
        bad_names = threshold_counts[threshold_counts != len(expected_thresholds)].index.astype(str).tolist()
        fail(f"seeded threshold cache does not contain ten thresholds for each run: {bad_names[:5]}")

    cached_seed_rows = (
        cached[["name", "seed_index", "seed_label", "seed"]]
        .drop_duplicates()
        .sort_values("name")
        .reset_index(drop=True)
    )
    manifest_seed_by_run = (
        run_manifest[["name", "seed_index", "seed_label", "seed"]]
        .drop_duplicates()
        .sort_values("name")
        .reset_index(drop=True)
    )
    max_dataframe_difference(cached_seed_rows, manifest_seed_by_run, label="seeded cache seed metadata")

    ok(f"seeded manifests match cached threshold table; runs={len(run_manifest)}, seeds={len(seed_manifest)}")
    return run_manifest


def validate_seeded_cache(namespace: dict[str, Any], repo_root: Path, tolerance: float, *, require_seeded: bool) -> None:
    seeded_cache_path = repo_root / "tables" / "seeded" / "threshold_long.csv"
    if not seeded_cache_path.exists():
        if require_seeded:
            fail(f"seeded threshold cache is missing: {seeded_cache_path}")
        skip("seeded threshold cache is absent")
        return

    cached = pd.read_csv(seeded_cache_path)
    run_manifest = validate_seeded_manifests(repo_root, cached)

    try:
        recomputed = namespace["compute_seeded_threshold_long_df"]()
    except FileNotFoundError as exc:
        if require_seeded:
            fail(str(exc))
        skip(f"seeded run root unavailable: {exc}")
        return

    sort_columns = ["display_name", "seed_index", "threshold"]
    cached = cached.sort_values(sort_columns).reset_index(drop=True)
    recomputed = recomputed.sort_values(sort_columns).reset_index(drop=True)
    for column in SEEDED_COMPARE_COLUMNS:
        diff = (cached[column].astype(float) - recomputed[column].astype(float)).abs()
        max_diff = float(diff.max()) if len(diff) else 0.0
        if max_diff > tolerance:
            fail(f"seeded cache differs from recomputation in {column}: {max_diff:.3g}")

    expected_rows = run_manifest["name"].nunique() * 10
    if len(cached) != expected_rows:
        fail(f"seeded threshold cache has {len(cached)} rows; expected {expected_rows}")
    ok(f"seeded threshold cache matches recomputation; rows={len(cached)}")


def validate_seeded_figure_inputs(namespace: dict[str, Any]) -> None:
    scores = namespace["compute_plot_3_10_seeded_auc_scores"](include_ground_truth=False)
    if scores.empty:
        fail("seeded AUC score table is empty")
    for column in NORMALISED_AUC_COLUMNS:
        values = pd.to_numeric(scores[column], errors="coerce")
        if values.isna().any():
            fail(f"{column} contains NaN values")
        if not bool(((values >= 0.0) & (values <= 1.0)).all()):
            fail(f"{column} has values outside [0, 1]")
        std_column = f"{column}_std"
        std_values = pd.to_numeric(scores[std_column], errors="coerce")
        if std_values.isna().any() or not bool((std_values >= 0.0).all()):
            fail(f"{std_column} contains invalid standard deviations")

    seeded_models = scores[~scores["display_name"].isin(["sup", "UDP0057"])]
    if not bool((seeded_models["seed_count"].astype(int) >= 2).all()):
        fail("one or more seeded model summaries has fewer than two seeds")

    all_scores = namespace["select_plot_3_10_all_models"](
        namespace["compute_plot_3_10_seeded_auc_scores"](),
    )
    for abundance_source in ("model", "ground_truth"):
        ecdf_df = namespace["build_plot_3_10_seeded_abundance_ecdf_df"](
            all_scores,
            abundance_source=abundance_source,
        )
        if ecdf_df.empty:
            fail(f"{abundance_source} ECDF table is empty")
        for display_name, group in ecdf_df.groupby("display_name", sort=False):
            abundance = pd.to_numeric(group["abundance"], errors="coerce").dropna()
            if abundance.empty:
                fail(f"{abundance_source} ECDF observations for {display_name} are empty")
            if not bool((abundance > 0).all()):
                fail(f"{abundance_source} ECDF observations for {display_name} have non-positive abundance")
            if "cumulative_probability" in group.columns:
                curve = group[["abundance", "cumulative_probability"]].dropna().sort_values("abundance")
                probabilities = curve["cumulative_probability"].to_numpy(dtype=float)
                if not bool(((probabilities >= 0.0) & (probabilities <= 1.0)).all()):
                    fail(f"{abundance_source} ECDF curve for {display_name} has probabilities outside [0, 1]")
                if np.any(np.diff(probabilities) < -1e-12):
                    fail(f"{abundance_source} ECDF curve for {display_name} is not monotonic")
                if abs(float(probabilities[-1]) - 1.0) > 1e-12:
                    fail(f"{abundance_source} ECDF curve for {display_name} does not end at 1")
    ok("seeded AUC and ECDF figure inputs satisfy numerical invariants")


def validate_tables(repo_root: Path) -> None:
    issues: list[str] = []
    table_paths = sorted((repo_root / "tables").rglob("*.csv"))
    if not table_paths:
        fail("no CSV tables found")
    optional_all_null_columns = {"seed_index", "seed_label"}
    for path in table_paths:
        df = pd.read_csv(path)
        rel = path.relative_to(repo_root)
        if df.empty:
            issues.append(f"{rel}: empty")
        unnamed = [column for column in df.columns if str(column).startswith("Unnamed")]
        if unnamed:
            issues.append(f"{rel}: unnamed columns {unnamed}")
        numeric = df.select_dtypes(include=[np.number])
        for column in numeric.columns:
            values = numeric[column].to_numpy(dtype=float)
            inf_mask = np.isinf(values)
            if not inf_mask.any():
                continue
            allowed_open_histogram_bin = (
                rel == Path("tables/weights/weight_delta_histogram.csv")
                and column == "bin_right"
                and "bin_index" in df.columns
                and bool((df.loc[inf_mask, "bin_index"] == df["bin_index"].max()).all())
            )
            if not allowed_open_histogram_bin:
                issues.append(f"{rel}: column {column} contains inf/-inf")
        all_null = [
            column
            for column in df.columns
            if df[column].isna().all() and column not in optional_all_null_columns
        ]
        if all_null:
            issues.append(f"{rel}: all-null columns {all_null}")
    if issues:
        fail("table integrity issues:\n" + "\n".join(issues))
    ok(f"table integrity checks passed for {len(table_paths)} CSV files")


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=None,
        help="Repository root. Defaults to searching upward from the current directory.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-9,
        help="Absolute tolerance for regenerated metric-table comparisons.",
    )
    parser.add_argument(
        "--skip-correlation-regeneration",
        action="store_true",
        help="Skip the slower end-to-end GenerateCorrelationCsv.sh check.",
    )
    parser.add_argument(
        "--skip-seeded",
        action="store_true",
        help="Skip seeded run-root/cache/figure-input checks.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_argument_parser().parse_args(argv)
    repo_root = args.repo_root.resolve() if args.repo_root else find_repo_root(Path.cwd())
    validate_notebook_syntax(repo_root)
    namespace = load_plot_namespace(repo_root)
    validate_threshold_semantics(namespace, repo_root)
    validate_notebook_concordance_matches_correlation(namespace, args.tolerance)
    if args.skip_correlation_regeneration:
        skip("end-to-end correlation regeneration")
    else:
        validate_generate_correlation(repo_root, args.tolerance)
    if args.skip_seeded:
        skip("seeded metric validations")
    else:
        validate_seeded_cache(namespace, repo_root, args.tolerance, require_seeded=True)
        validate_seeded_figure_inputs(namespace)
    validate_tables(repo_root)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationFailure as exc:
        print(f"[fail] {exc}", file=sys.stderr)
        raise SystemExit(1)
