#!/usr/bin/env python3
"""Reproduce spatial convergence, time-step sensitivity and solver comparisons.

Run in the deal.II environment after building. No shell commands are constructed
from user input. Only named output files are replaced; folders are never deleted.
"""

import argparse
import base64
import csv
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import shutil
import xml.etree.ElementTree as ET
import zlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_table(path, columns, data):
    with path.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(data)
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    lines += ["| " + " | ".join(str(row[c]) for c in columns) + " |" for row in data]
    path.with_suffix(".md").write_text("\n".join(lines) + "\n")
    selections = {
        "convergence": ["case", "N", "L2", "H1", "p_L2", "p_H1"],
        "time-sensitivity": ["case", "N", "L2_change", "H1_change"],
        "energy-summary": ["case", "initial", "final", "minimum", "max_relative_change"],
        "solver-comparison": ["solver", "preconditioner", "converged", "mean_iterations", "max_residual", "solve_ms", "setup_ms"],
    }
    if path.name in selections:
        selected = selections[path.name]
        labels = {"p_L2": "$p_{L^2}$", "p_H1": "$p_{H^1}$", "L2": "$L^2$", "H1": "$H^1$",
                  "L2_change": "$L^2$ change", "H1_change": "$H^1$ change",
                  "max_relative_change": "max. rel. change", "preconditioner": "prec.",
                  "mean_iterations": "mean it.", "max_residual": "max. residual",
                  "solve_ms": "solve ms", "setup_ms": "setup ms"}
        latex = [r"\begin{center}\small", r"\begin{tabular}{" + "l" * len(selected) + "}",
                 r"\toprule", " & ".join(labels.get(c, c) for c in selected) + r"\\", r"\midrule"]
        latex += [" & ".join(str(row[c]) for c in selected) + r"\\" for row in data]
        latex += [r"\bottomrule", r"\end{tabular}\end{center}"]
        path.with_suffix(".tex").write_text("\n".join(latex) + "\n")


def vtk_array(element):
    # deal.II writes separately base64-encoded zlib header and compressed data.
    encoded = "".join(element.text.split())
    nblocks = int(np.frombuffer(base64.b64decode(encoded[:8])[:4], dtype="<u4")[0])
    header_chars = 4 * math.ceil((3 + nblocks) * 4 / 3)
    header = np.frombuffer(base64.b64decode(encoded[:header_chars]), dtype="<u4")
    compressed = base64.b64decode(encoded[header_chars:])
    offset = 0
    blocks = []
    for size in header[3:]:
        blocks.append(zlib.decompress(compressed[offset:offset + size]))
        offset += int(size)
    types = {"Float32": "<f4", "Float64": "<f8", "UInt8": "u1", "Int32": "<i4"}
    return np.frombuffer(b"".join(blocks), dtype=types[element.attrib["type"]])


def check_output(directory, final_time):
    data = rows(directory / "diagnostics.csv")
    assert math.isclose(float(data[-1]["time"]), final_time, abs_tol=1.e-12)
    assert all(math.isfinite(float(row[key])) for row in data for key in
               ("l2_error", "h1_seminorm", "energy"))
    entries = ET.parse(directory / "wave_solution.pvd").findall(".//DataSet")
    assert math.isclose(float(entries[-1].attrib["timestep"]), final_time, abs_tol=1.e-12)
    for entry in entries:
        tree = ET.parse(directory / entry.attrib["file"])
        assert set(vtk_array(tree.find(".//Cells/DataArray[@Name='types']"))) == {9}
        points = vtk_array(tree.find(".//Points/DataArray")).reshape(-1, 3)
        boundary = np.any(np.isclose(points[:, :2], 0) | np.isclose(points[:, :2], 1), axis=1)
        for field in ("u", "velocity", "acceleration"):
            values = vtk_array(tree.find(f".//PointData/DataArray[@Name='{field}']"))
            assert np.all(np.isfinite(values)) and np.max(np.abs(values[boundary])) < 1.e-12
    assert all(r["converged"] == "1" for r in rows(directory / "solver_stats.csv"))
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path("output/validation/review"))
    parser.add_argument("--reuse", action="store_true", help="Regenerate tables/figures from existing runs")
    args = parser.parse_args()
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    assets = Path("submission/assets")
    tables = assets / "tables"
    figures = assets / "figures"
    tables.mkdir(exist_ok=True)
    env = dict(os.environ, WAVE_CHECK_MATRIX="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

    def run(name, case, n, dt, final=1, solver="cg", prec="ssor", degree=1):
        directory = root / name
        directory.mkdir(parents=True, exist_ok=True)
        command = [str(Path("build/wave_fem").resolve()), str(n), str(degree), str(final),
                   str(dt), str(round(final / dt)), str(directory), case, solver, prec]
        if not args.reuse:
            print("Running", name, flush=True)
            with (directory / "run.log").open("w") as log:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                return directory, None
        metadata = rows(directory / "parameters.csv")[0]
        assert (metadata["case"], int(metadata["N"]), int(metadata["degree"]),
                float(metadata["dt"]), float(metadata["T"]), metadata["solver"], metadata["preconditioner"]) == (
                case, n, degree, dt, final, solver, prec)
        return directory, check_output(directory, final)

    convergence = []
    sensitivity = []
    energies = []
    checks = []
    for case in ("free", "forced"):
        previous = None
        case_results = []
        for n in (12, 24, 48):
            directory, data = run(f"{case}/space_N{n}", case, n, 0.00125)
            assert data is not None, directory
            _, finer = run(f"{case}/half_dt_N{n}", case, n, 0.000625)
            assert finer is not None
            final = data[-1]
            h = math.sqrt(2) / n
            errors = [float(final[k]) for k in ("l2_error", "h1_seminorm")]
            eoc = [math.log(previous[1][i] / errors[i]) / math.log(previous[0] / h)
                   if previous else float("nan") for i in range(2)]
            record = dict(case=case, N=n, h=f"{h:.8g}", dt=0.00125,
                          L2=f"{errors[0]:.8e}", H1=f"{errors[1]:.8e}",
                          p_L2=f"{eoc[0]:.4f}" if previous else "-",
                          p_H1=f"{eoc[1]:.4f}" if previous else "-")
            convergence.append(record)
            case_results.append(record)
            previous = h, errors
            # Norm differences are a sensitivity indicator, not a rigorous error bound.
            changes = [abs(float(finer[-1][key]) - errors[i]) / errors[i]
                       for i, key in enumerate(("l2_error", "h1_seminorm"))]
            sensitivity.append(dict(case=case, N=n, L2_change=f"{changes[0]:.6e}",
                                    H1_change=f"{changes[1]:.6e}"))
            assert max(changes) < 0.02, "Time-step contamination exceeds 2%"
            if n == 48:
                initial = float(data[0]["energy"])
                energy = [float(r["energy"]) for r in data]
                drift = max(abs(e - initial) / initial for e in energy)
                energies.append(dict(case=case, initial=f"{initial:.12g}",
                    final=f"{energy[-1]:.12g}", minimum=f"{min(energy):.12g}",
                    maximum=f"{max(energy):.12g}", max_relative_change=f"{drift:.6e}"))
                if case == "free":
                    assert drift < 1.e-8
                subprocess.run([os.sys.executable, "scripts/plot_diagnostics.py",
                                str(directory / "diagnostics.csv")], check=True)
                for kind in ("errors", "energy", "energy_change"):
                    shutil.copyfile(directory / f"{kind}.png", figures / f"{case}-time-{kind}.png")
                for row in rows(directory / "matrix_checks.csv"):
                    checks.append(dict(case=case, **row))
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        hs = np.array([float(r["h"]) for r in case_results])
        for key, slope, color in (("L2", 2, "tab:blue"), ("H1", 1, "tab:orange")):
            es = np.array([float(r[key]) for r in case_results])
            ax.loglog(hs, es, "o-", color=color, label=f"absolute {key} {'seminorm' if key == 'H1' else 'norm'}")
            ax.loglog(hs, es[-1] * (hs / hs[-1]) ** slope, "--", color=color, label=f"reference slope {slope}")
        ax.set(xlabel="h = maximum cell diameter = sqrt(2)/N", ylabel="absolute error at T = 1",
               title=f"{case.capitalize()}: Q1 spatial convergence, dt = 0.00125")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(figures / f"{case}-convergence.png", dpi=200)
        plt.close(fig)

    benchmarks = []
    benchmark_diagnostics = []
    for solver in ("cg", "gmres"):
        for prec in ("ssor", "jacobi"):
            trials = []
            for trial in range(3):
                directory, data = run(f"solvers/{solver}-{prec}-{trial}", "forced", 48, 0.005,
                                      final=0.1, solver=solver, prec=prec)
                stats = rows(directory / "solver_stats.csv")
                steps = [r for r in stats if r["system"] == "Newmark acceleration"]
                success = data is not None and len(steps) == 20
                if success:
                    benchmark_diagnostics.append([float(data[-1][k]) for k in ("l2_error", "h1_seminorm", "energy")])
                trials.append(dict(success=success, steps=steps))
            steps = [row for trial in trials for row in trial["steps"]]
            benchmarks.append(dict(solver=solver, preconditioner=prec,
                converged=f"{sum(t['success'] for t in trials)}/3",
                mean_iterations=f"{statistics.mean(int(r['iterations']) for r in steps):.2f}" if steps else "-",
                max_iterations=max((int(r["iterations"]) for r in steps), default=0),
                max_residual=f"{max((float(r['residual']) for r in steps), default=math.inf):.6e}",
                max_relative_residual=f"{max((float(r['relative_residual']) for r in steps), default=math.inf):.6e}",
                solve_ms=f"{1000 * statistics.median(sum(float(r['solve_seconds']) for r in t['steps']) for t in trials):.5f}",
                setup_ms=f"{1000 * statistics.median(sum(float(r['setup_seconds']) for r in t['steps']) for t in trials):.5f}"))
    assert all(row["converged"] == "3/3" for row in benchmarks)
    assert np.allclose(benchmark_diagnostics, benchmark_diagnostics[0], rtol=1.e-7, atol=1.e-10)
    _, higher_degree = run("smoke_Q2", "forced", 8, 0.005, final=0.05, degree=2)
    assert higher_degree is not None
    executable = str(Path("build/wave_fem").resolve())
    valid = ["8", "1", "0.05", "0.005", "10", str(root / "invalid"), "free", "cg", "ssor"]
    for index, value in ((0, "0"), (0, "-1"), (0, "3abc"), (1, "0"),
                         (2, "nan"), (3, "0"), (3, "0.03"), (4, "0"),
                         (6, "unknown"), (7, "unknown"), (8, "unknown")):
        arguments = valid.copy()
        arguments[index] = value
        result = subprocess.run([executable, *arguments], env=env, capture_output=True, timeout=10)
        assert result.returncode == 1 and (b"Error:" in result.stderr or
               b"Unknown validation_case" in result.stderr), (arguments, result.stderr)
    for name, data in (("convergence", convergence), ("time-sensitivity", sensitivity),
                       ("energy-summary", energies), ("solver-comparison", benchmarks), ("matrix-checks", checks)):
        write_table(tables / name, list(data[0]), data)
    environment_file = root / "benchmark-environment.txt"
    if not args.reuse:
        environment_file.write_text(
            f"Platform: {platform.platform()}\nProcessor: {platform.processor()}\n"
            "Reference toolchain: course deal.II 9.5.1, CMake Release.\n"
            "Execution: serial, OMP/OPENBLAS threads=1.\n"
            "Benchmark: forced Q1 N=48 (2401 DoFs), T=0.1, dt=0.005; three repetitions.\n"
            "Timings: median total for 20 Newmark solves; setup separate; initial mass solve excluded.\n"
            "Residual: maximum recomputed Euclidean norm over Newmark solves/repetitions.\n"
            "Stats for every solve, including the initial mass solve, are in each solver_stats.csv.\n")
    shutil.copyfile(environment_file, tables / environment_file.name)
    print("All numerical, matrix, boundary, VTK and comparison checks passed.", flush=True)


if __name__ == "__main__":
    main()
