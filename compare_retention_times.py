#!/usr/bin/env python3
"""Compare apex retention times for matching theoretical m/z values.

For each theoretical m/z, the script calculates Q1, Q3, median, and mean
retention-time summaries across all input CSV files. Optionally, each value
can be labelled as lower_outlier (< Q1), within_iqr (Q1 through Q3),
upper_outlier (> Q3), or missing. The output is transposed: theoretical m/z
values are columns and measurements are rows.

Set ``CSV_FOLDER`` near the top of this file, then run the script without
command-line arguments. Only the Python standard library is required.
"""

from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path
from statistics import fmean, median
from typing import Sequence


# ---------------------------------------------------------------------------
# USER SETTINGS: change CSV_FOLDER to the folder containing the CSV files.
# Examples:
#   Windows: CSV_FOLDER = Path(r"C:\data\retention_times")
#   macOS/Linux: CSV_FOLDER = Path("/home/user/data/retention_times")
# The default below uses the folder containing this Python script.
# ---------------------------------------------------------------------------
CSV_FOLDER = Path("/home/haotian/IsotopePairFinder/mzML/Mass_Accuracy")
OUTPUT_FILENAME = "retention_time_comparison.csv"
INCLUDE_OUTLIER_LABELS = False

MZ_COLUMN = "theoretical_mz"
RT_COLUMN = "apex_retention_time_min"


def percentile(values: Sequence[float], probability: float) -> float:
    """Return a linearly interpolated percentile (the common type-7 method)."""
    if not values:
        raise ValueError("Cannot calculate a percentile of an empty sequence")
    if not 0.0 <= probability <= 1.0:
        raise ValueError("Percentile probability must be between 0 and 1")

    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower_index = math.floor(position)
    upper_index = math.ceil(position)
    if lower_index == upper_index:
        return ordered[lower_index]

    fraction = position - lower_index
    return ordered[lower_index] + (
        ordered[upper_index] - ordered[lower_index]
    ) * fraction


def classify_retention_time(value: float | None, q1: float, q3: float) -> str:
    if value is None:
        return "missing"
    if value < q1:
        return "lower_outlier"
    if value > q3:
        return "upper_outlier"
    return "within_iqr"


def parse_mz(raw_value: str, path: Path, line_number: int) -> Decimal:
    try:
        value = Decimal(raw_value.strip())
    except InvalidOperation as exc:
        raise ValueError(
            f"{path}: line {line_number}: invalid {MZ_COLUMN} value "
            f"{raw_value!r}"
        ) from exc
    if not value.is_finite():
        raise ValueError(
            f"{path}: line {line_number}: {MZ_COLUMN} must be finite"
        )
    return value


def parse_retention_time(
    raw_value: str, path: Path, line_number: int
) -> float | None:
    raw_value = raw_value.strip()
    if not raw_value:
        return None
    try:
        value = float(raw_value)
    except ValueError as exc:
        raise ValueError(
            f"{path}: line {line_number}: invalid {RT_COLUMN} value "
            f"{raw_value!r}"
        ) from exc
    if not math.isfinite(value):
        raise ValueError(
            f"{path}: line {line_number}: {RT_COLUMN} must be finite"
        )
    return value


def read_input_csv(path: Path) -> tuple[dict[Decimal, float | None], dict[Decimal, str]]:
    """Read one input file and return RT values plus original m/z spellings."""
    values: dict[Decimal, float | None] = {}
    mz_display: dict[Decimal, str] = {}

    with path.open("r", encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: the CSV file has no header")

        missing_columns = [
            name for name in (MZ_COLUMN, RT_COLUMN) if name not in reader.fieldnames
        ]
        if missing_columns:
            raise ValueError(
                f"{path}: missing required column(s): {', '.join(missing_columns)}"
            )

        for line_number, row in enumerate(reader, start=2):
            raw_mz = (row.get(MZ_COLUMN) or "").strip()
            if not raw_mz:
                continue

            mz = parse_mz(raw_mz, path, line_number)
            if mz in values:
                raise ValueError(
                    f"{path}: line {line_number}: duplicate {MZ_COLUMN} {raw_mz!r}"
                )

            values[mz] = parse_retention_time(
                row.get(RT_COLUMN) or "", path, line_number
            )
            mz_display[mz] = raw_mz

    return values, mz_display


def has_required_columns(path: Path) -> bool:
    """Return whether a CSV has the two columns used by this comparison."""
    with path.open("r", encoding="utf-8-sig", newline="") as input_file:
        reader = csv.reader(input_file)
        try:
            fieldnames = next(reader)
        except StopIteration:
            return False
    return MZ_COLUMN in fieldnames and RT_COLUMN in fieldnames


def find_input_csvs(folder: Path, output_path: Path) -> tuple[list[Path], list[Path]]:
    """Find CSVs directly inside folder, excluding output and unrelated CSVs."""
    output_resolved = output_path.resolve()
    candidates = sorted(
        (
            path
            for path in folder.iterdir()
            if path.is_file()
            and path.suffix.lower() == ".csv"
            and path.resolve() != output_resolved
        ),
        key=lambda path: path.name.casefold(),
    )

    input_paths: list[Path] = []
    skipped_paths: list[Path] = []
    for path in candidates:
        if has_required_columns(path):
            input_paths.append(path)
        else:
            skipped_paths.append(path)
    return input_paths, skipped_paths


def unique_file_labels(paths: Sequence[Path]) -> list[str]:
    """Return readable, unique labels for output column names."""
    counts: Counter[str] = Counter()
    labels: list[str] = []
    for path in paths:
        counts[path.name] += 1
        occurrence = counts[path.name]
        labels.append(path.name if occurrence == 1 else f"{path.name}_{occurrence}")
    return labels


def format_number(value: float | None) -> str:
    return "" if value is None else format(value, ".12g")


def compare_files(
    input_paths: Sequence[Path],
    output_path: Path,
    include_outlier_labels: bool = INCLUDE_OUTLIER_LABELS,
) -> tuple[int, Counter[str]]:
    file_values: list[dict[Decimal, float | None]] = []
    mz_display: dict[Decimal, str] = {}

    for path in input_paths:
        values, displays = read_input_csv(path)
        file_values.append(values)
        for mz, display in displays.items():
            mz_display.setdefault(mz, display)

    labels = unique_file_labels(input_paths)
    label_counts: Counter[str] = Counter()
    comparison_rows: list[dict[str, str]] = []
    sorted_mz = sorted(mz_display)

    for mz in sorted_mz:
        retention_times = [values.get(mz) for values in file_values]
        observed = [value for value in retention_times if value is not None]

        row: dict[str, str] = {MZ_COLUMN: mz_display[mz]}
        if observed:
            q1 = percentile(observed, 0.25)
            q3 = percentile(observed, 0.75)
            row["q1_retention_time_min"] = format_number(q1)
            row["q3_retention_time_min"] = format_number(q3)
            row["median_retention_time_min"] = format_number(median(observed))
            row["mean_retention_time_min"] = format_number(fmean(observed))
        else:
            q1 = q3 = math.nan
            row["q1_retention_time_min"] = ""
            row["q3_retention_time_min"] = ""
            row["median_retention_time_min"] = ""
            row["mean_retention_time_min"] = ""

        for label, value in zip(labels, retention_times):
            if observed:
                classification = classify_retention_time(value, q1, q3)
            else:
                classification = "missing"
            row[f"{label}__{RT_COLUMN}"] = format_number(value)
            row[f"{label}__outlier_label"] = classification
            label_counts[classification] += 1

        comparison_rows.append(row)

    with output_path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file)
        mz_headers = [mz_display[mz] for mz in sorted_mz]
        writer.writerow(["source", "measurement", *mz_headers])
        for fieldname in (
            "q1_retention_time_min",
            "q3_retention_time_min",
            "median_retention_time_min",
            "mean_retention_time_min",
        ):
            writer.writerow(
                ["summary", fieldname, *(row[fieldname] for row in comparison_rows)]
            )
        for label in labels:
            writer.writerow(
                [
                    label,
                    RT_COLUMN,
                    *(row[f"{label}__{RT_COLUMN}"] for row in comparison_rows),
                ]
            )
            if include_outlier_labels:
                writer.writerow(
                    [
                        label,
                        "outlier_label",
                        *(row[f"{label}__outlier_label"] for row in comparison_rows),
                    ]
                )

    return len(sorted_mz), label_counts


def main() -> int:
    folder = CSV_FOLDER.expanduser()
    output_path = folder / OUTPUT_FILENAME

    if not folder.is_dir():
        print(
            f"error: CSV_FOLDER does not exist or is not a folder: {folder}",
            file=sys.stderr,
        )
        return 1

    try:
        input_paths, skipped_paths = find_input_csvs(folder, output_path)
        if len(input_paths) < 2:
            raise ValueError(
                f"found {len(input_paths)} compatible CSV file(s) in {folder}; "
                "at least 2 are required"
            )

        print(f"Comparing {len(input_paths)} CSV files from {folder}")
        for path in skipped_paths:
            print(f"Skipping {path.name}: required columns not found")

        mz_count, label_counts = compare_files(
            input_paths,
            output_path,
            include_outlier_labels=INCLUDE_OUTLIER_LABELS,
        )
    except (OSError, ValueError, csv.Error) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {mz_count} theoretical m/z rows to {output_path}")
    if INCLUDE_OUTLIER_LABELS:
        print(
            "Labels: "
            + ", ".join(
                f"{name}={label_counts[name]}"
                for name in (
                    "lower_outlier",
                    "within_iqr",
                    "upper_outlier",
                    "missing",
                )
            )
        )
    else:
        print("Outlier label rows are disabled")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
