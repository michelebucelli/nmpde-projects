#!/usr/bin/env python3
"""Create a compact Markdown summary from WaveFEM diagnostics files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_final_row(csv_file: Path) -> dict[str, str]:
    with csv_file.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ValueError(f"{csv_file} has no data rows")

    return rows[-1]


def infer_case_and_run(csv_file: Path, root: Path) -> tuple[str, str]:
    relative = csv_file.parent.relative_to(root)
    parts = relative.parts

    if len(parts) >= 2:
        return parts[0], "/".join(parts[1:])

    return "default", str(relative)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "root",
        nargs="?",
        default="output/validation/suite",
        help="Folder containing validation runs.",
    )
    parser.add_argument(
        "--output",
        default="submission/assets/tables/validation-summary.md",
        help="Markdown file to write.",
    )
    args = parser.parse_args()

    root = Path(args.root)
    output = Path(args.output)

    csv_files = sorted(root.glob("**/diagnostics.csv"))
    if not csv_files:
        raise SystemExit(f"No diagnostics.csv files found under {root}")

    lines = [
        "# Validation Summary",
        "",
        "| case | run | final time | L2 error | H1 seminorm | energy |",
        "|---|---:|---:|---:|---:|---:|",
    ]

    for csv_file in csv_files:
        row = read_final_row(csv_file)
        case, run_name = infer_case_and_run(csv_file, root)
        lines.append(
            "| "
            f"{case} | "
            f"{run_name} | "
            f"{float(row['time']):.6g} | "
            f"{float(row['l2_error']):.6e} | "
            f"{float(row['h1_seminorm']):.6e} | "
            f"{float(row['energy']):.6e} |"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
