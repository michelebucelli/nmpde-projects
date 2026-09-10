#!/usr/bin/env python3
"""Measure the first centre-probe peak for several Newmark time steps."""

import csv
import math
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "output/validation/dispersion"
TABLE = ROOT / "submission/assets/tables/dispersion-check.csv"
MARKDOWN = ROOT / "submission/assets/tables/dispersion-check.md"
FIGURE = ROOT / "submission/assets/figures/dispersion-check.png"
LATEX = ROOT / "submission/assets/tables/dispersion-check.tex"


def read_probe(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return (np.array([float(row["time"]) for row in rows]),
            np.array([float(row["numerical"]) for row in rows]),
            np.array([float(row["exact"]) for row in rows]))


def peak_time(times, values, expected):
    window = np.flatnonzero((times > 0.75 * expected) & (times < 1.25 * expected))
    index = window[np.argmax(values[window])]
    if index == 0 or index == len(times) - 1:
        return times[index]
    # A quadratic through three samples removes most of the output-grid error.
    coefficients = np.polyfit(times[index - 1:index + 2],
                              values[index - 1:index + 2], 2)
    return float(-coefficients[1] / (2.0 * coefficients[0]))


def main():
    RUN_ROOT.mkdir(exist_ok=True)
    executable = ROOT / "build/wave_fem"
    omega = math.sqrt(2.0) * math.pi
    exact_period = 2.0 * math.pi / omega
    records = []
    curves = []
    for dt in (0.02, 0.01, 0.005):
        directory = RUN_ROOT / f"dt_{dt:g}"
        directory.mkdir(parents=True, exist_ok=True)
        command = [str(executable), "192", "1", "2.0", str(dt),
                   str(round(2.0 / dt)), str(directory), "free"]
        subprocess.run(command, cwd=ROOT, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        times, numerical, exact = read_probe(directory / "probe.csv")
        first_peak = peak_time(times, numerical, exact_period)
        lag = first_peak - exact_period
        records.append({"dt": f"{dt:.5f}", "exact_period": f"{exact_period:.8f}",
                        "numerical_peak": f"{first_peak:.8f}",
                        "peak_lag": f"{lag:.8e}",
                        "phase_lag_rad": f"{omega * lag:.8e}"})
        if dt == 0.01:
            curves = (times, numerical, exact)

    TABLE.parent.mkdir(parents=True, exist_ok=True)
    columns = ["dt", "exact_period", "numerical_peak", "peak_lag", "phase_lag_rad"]
    with TABLE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader(); writer.writerows(records)
    lines = ["| dt | exact period | numerical first peak | peak lag | phase lag (rad) |",
             "|---:|---:|---:|---:|---:|"]
    lines += [f"| {r['dt']} | {r['exact_period']} | {r['numerical_peak']} | "
              f"{r['peak_lag']} | {r['phase_lag_rad']} |" for r in records]
    MARKDOWN.write_text("\n".join(lines) + "\n")
    latex_lines = [
        r"\begin{center}",
        r"\small",
        r"\begin{tabular}{rrrrr}",
        r"\toprule",
        r"$\Delta t$ & exact period & numerical peak & peak lag & phase lag (rad) \\",
        r"\midrule",
    ]
    for record in records:
        latex_lines.append(
            f"{record['dt']} & {record['exact_period']} & "
            f"{record['numerical_peak']} & {record['peak_lag']} & "
            f"{record['phase_lag_rad']} \\\\")
    latex_lines += [r"\bottomrule", r"\end{tabular}", r"\end{center}"]
    LATEX.write_text("\n".join(latex_lines) + "\n")

    times, numerical, exact = curves
    mask = times <= 1.2 * exact_period
    plt.figure(figsize=(6.4, 3.8), dpi=180)
    plt.plot(times[mask], exact[mask], label="exact", color="#102C53")
    plt.plot(times[mask], numerical[mask], "--", label="numerical, dt=0.01",
             color="#4DC9FF")
    plt.axvline(exact_period, color="#102C53", alpha=0.3, linewidth=1)
    plt.xlabel("time")
    plt.ylabel("u(0.5, 0.5, t)")
    plt.title("Empirical phase check at the centre probe")
    plt.grid(alpha=0.25); plt.legend(); plt.tight_layout()
    plt.savefig(FIGURE); plt.close()
    print(f"Wrote {TABLE}, {LATEX} and {FIGURE}")


if __name__ == "__main__":
    main()
