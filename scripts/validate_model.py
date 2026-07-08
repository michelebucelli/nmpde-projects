#!/usr/bin/env python3
"""Model validator for this project.

Checks the theta-method mass-balance identity: for homogeneous Neumann BCs,
taking v_h = 1 (which lies in the FE space) in the fully discrete scheme
yields exactly

    integral_Omega u_h(T) = sum_n dt_n [theta * integral f^(n+1) + (1-theta) * integral f^n]

`Heat::run()` writes both sides of this identity to a two-column CSV
(`computed_mass,expected_mass`) after the run completes; this script compares
them and fails if the relative error exceeds the given tolerance.
"""

import argparse
import csv
import sys
from pathlib import Path


def read_masses(path: Path):
    with path.open(newline="") as fh:
        reader = csv.reader(fh)
        next(reader)
        row = next(reader)
    return float(row[0]), float(row[1])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", "-r", required=True, help="mass_balance.csv produced by the solver")
    p.add_argument("--tolerance", "-t", type=float, default=1e-6)
    args = p.parse_args()

    path = Path(args.results)
    if not path.exists():
        print(f"Results file not found: {path}", file=sys.stderr)
        sys.exit(2)

    computed_mass, expected_mass = read_masses(path)
    rel_err = abs(computed_mass - expected_mass) / max(abs(expected_mass), 1e-12)

    print(f"Computed mass:  {computed_mass:.12g}")
    print(f"Expected mass:  {expected_mass:.12g}  (theta-method mass balance)")
    print(f"Relative error: {rel_err:.6g} (tolerance {args.tolerance})")

    if rel_err > args.tolerance:
        print("Validation FAILED", file=sys.stderr)
        sys.exit(5)

    print("Validation PASSED")
    sys.exit(0)


if __name__ == "__main__":
    main()
