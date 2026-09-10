#!/usr/bin/env python3
"""Plot WaveFEM diagnostic CSV files.

The solver writes output/diagnostics.csv with time, error and energy. This
script creates quick PNG plots for the report/presentation. It is intentionally
small and optional: the C++ solver does not depend on Python.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_csv(path: Path) -> dict[str, list[float]]:
    data: dict[str, list[float]] = {
        "time": [],
        "l2_error": [],
        "h1_seminorm": [],
        "energy": [],
    }

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            for key in data:
                data[key].append(float(row[key]))

    return data


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csv_file",
        nargs="?",
        default="output/diagnostics.csv",
        help="Path to diagnostics.csv produced by wave_fem.",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv_file)
    output_dir = csv_path.parent

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required for plots. Install with:\n"
            "python -m pip install --user matplotlib"
        ) from exc

    data = read_csv(csv_path)

    plt.figure()
    plt.semilogy(data["time"], data["l2_error"], label="L2 error")
    plt.semilogy(data["time"], data["h1_seminorm"], label="H1 seminorm error")
    plt.xlabel("time")
    plt.ylabel("absolute spatial error norm at time t")
    plt.title("Displacement error evolution (fixed mesh and time step)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "errors.png", dpi=200)
    plt.close()

    plt.figure()
    initial_energy = data["energy"][0]
    plt.plot(data["time"], data["energy"], label="computed energy")
    plt.axhline(initial_energy, color="black", linestyle="--", label="initial energy")
    plt.xlabel("time")
    plt.ylabel("discrete energy")
    plt.ylim(0, max(1.e-15, 1.1 * max(data["energy"])))
    plt.ticklabel_format(axis="y", style="plain", useOffset=False)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_dir / "energy.png", dpi=200)
    plt.close()

    plt.figure()
    scale = initial_energy if initial_energy != 0 else 1.0
    plt.plot(data["time"], [(e - initial_energy) / scale for e in data["energy"]])
    plt.xlabel("time")
    plt.ylabel("(E(t) - E(0)) / E(0)" if initial_energy != 0 else "E(t) - E(0)")
    plt.title("Relative energy change (drift only for the free case)"
              if initial_energy != 0 else "Absolute energy change: zero initial energy")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_dir / "energy_change.png", dpi=200)
    plt.close()

    print(f"Wrote {output_dir / 'errors.png'}")
    print(f"Wrote {output_dir / 'energy.png'}")
    print(f"Wrote {output_dir / 'energy_change.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
