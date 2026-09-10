"""Derive the small, auditable tables consumed by report.tex from benchmark CSV."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "result" / "benchmark_results.csv"
OUT = ROOT / "report" / "generated"
METHODS = ["none", "jacobi", "ssor", "ilu", "amg"]


def load() -> list[dict[str, str]]:
    with SOURCE.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def find(rows: list[dict[str, str]], np: int, p: int, r: int, method: str) -> dict[str, str]:
    return next(
        row for row in rows
        if int(row["mpi_procs"]) == np
        and int(float(row["p"])) == p
        and int(row["refinements"]) == r
        and row["preconditioner"] == method
    )


def write(name: str, columns: list[str], records: list[list[object]]) -> None:
    path = OUT / name
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(records)


def main() -> None:
    rows = load()
    OUT.mkdir(parents=True, exist_ok=True)

    write("iterations_by_p.dat", ["p", *METHODS], [
        [p, *[find(rows, 1, p, 5, method)["iterations"] for method in METHODS]]
        for p in range(1, 6)
    ])
    write("condition_by_p.dat", ["p", *METHODS], [
        [p, *[find(rows, 1, p, 5, method)["condition_number"] for method in METHODS]]
        for p in range(1, 6)
    ])
    write("iterations_by_ref.dat", ["r", *METHODS], [
        [r, *[find(rows, 1, 5, r, method)["iterations"] for method in METHODS]]
        for r in range(3, 6)
    ])

    serial = {method: float(find(rows, 1, 5, 5, method)["total_s"]) for method in METHODS}
    write("speedup_by_np.dat", ["np", *METHODS], [
        [np, *[serial[method] / float(find(rows, np, 5, 5, method)["total_s"])
              for method in METHODS]]
        for np in range(1, 5)
    ])


if __name__ == "__main__":
    main()
