from __future__ import annotations

import csv
import math
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULT = ROOT / "result"
PLOTS = RESULT / "plots"
EXECUTABLE = ROOT / "build" / "elliptic"
SETUP = ROOT / "setup_modules.sh"
CSV_FILE = RESULT / "benchmark_results.csv"

MPI_PROCS = [1, 2, 3, 4]
PRECONDITIONERS = ["none", "jacobi", "ssor", "ilu", "amg"]

RESULT_RE = re.compile(
    r"^result\s+(\S+)\s+(\d+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)(?:\s+([0-9.eE+-]+))?"
)


def run_shell(command: str, log_name: str, cwd: Path = ROOT) -> int:
    completed = subprocess.run(
        ["bash", "-lc", command],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (RESULT / log_name).write_text(completed.stdout, encoding="utf-8")
    return completed.returncode


def shell_output(command: str) -> str:
    completed = subprocess.run(
        ["bash", "-lc", command],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    return completed.stdout.strip()


def build_code() -> int:
    return run_shell(f"source {SETUP} && cmake --build build", "build.log")


def mpi_launcher() -> str:
    return shell_output(
        f"source {SETUP} && "
        "for cmd in mpirun mpiexec srun; do command -v $cmd && break; done"
    )


def run_experiment(np: int, launcher: str) -> int:
    if np == 1:
        command = f"source {SETUP} && {EXECUTABLE}"
    elif launcher:
        flag = "-n" if Path(launcher).name == "srun" else "-np"
        command = f"source {SETUP} && {launcher} {flag} {np} {EXECUTABLE}"
    else:
        (RESULT / f"run_np{np}.log").write_text(
            "mpi launcher missing\n",
            encoding="utf-8",
        )
        return 127

    return run_shell(command, f"run_np{np}.log", RESULT)


def parse_run_log(np: int, return_code: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    current = {
        "status": "ok",
        "returncode": str(return_code),
        "mpi_procs": str(np),
        "p": "",
        "refinements": "",
        "balls": "",
        "cells": "",
        "dofs": "",
        "case_setup_s": "",
        "assembly_s": "",
    }

    for line in (RESULT / f"run_np{np}.log").read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith("mpi_procs "):
            current["mpi_procs"] = line.split()[1]
        elif line.startswith("p "):
            current["p"] = line.split()[1]
        elif line.startswith("refinements "):
            current["refinements"] = line.split()[1]
        elif line.startswith("balls "):
            current["balls"] = line.split()[1]
        elif line.startswith("cells "):
            current["cells"] = line.split()[1]
        elif line.startswith("dofs "):
            current["dofs"] = line.split()[1]
        elif line.startswith("case_setup "):
            current["case_setup_s"] = line.split()[1]
        elif line.startswith("assembly "):
            current["assembly_s"] = line.split()[1]
        elif match := RESULT_RE.match(line):
            rows.append(
                current
                | {
                    "preconditioner": match.group(1),
                    "iterations":     match.group(2),
                    "pc_setup_s":     match.group(3),
                    "solve_s":        match.group(4),
                    "total_s":        match.group(5),
                    "condition_number": match.group(6) or "",
                }
            )

    return rows


def write_csv(rows: list[dict[str, str]]) -> None:
    columns = [
        "status",
        "returncode",
        "mpi_procs",
        "p",
        "refinements",
        "balls",
        "cells",
        "dofs",
        "case_setup_s",
        "assembly_s",
        "preconditioner",
        "iterations",
        "pc_setup_s",
        "solve_s",
        "total_s",
        "condition_number"
    ]
    with CSV_FILE.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def row_for(rows: list[dict[str, str]], np: str, p: str, ref: str, prec: str):
    for row in rows:
        if (
            row["mpi_procs"] == np
            and row["p"] == p
            and row["refinements"] == ref
            and row["preconditioner"] == prec
        ):
            return row
    return None


def sorted_values(rows: list[dict[str, str]], key: str) -> list[str]:
    values = {row[key] for row in rows if row.get(key, "")}
    return sorted(values, key=float)


def first_value(rows: list[dict[str, str]], key: str) -> str:
    values = sorted_values(rows, key)
    return values[0] if values else ""


def last_value(rows: list[dict[str, str]], key: str) -> str:
    values = sorted_values(rows, key)
    return values[-1] if values else ""


def number_label(value: str | float) -> str:
    value = float(value)
    return str(int(value)) if value.is_integer() else f"{value:g}"


def positive_max(values: list[float], fallback: float) -> float:
    values = [value for value in values if value > 0]
    return max(values) if values else fallback


def svg_start(title: str, width: int = 860, height: int = 520) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>text{font-family:Arial,sans-serif;font-size:13px;fill:#1f2933}.title{font-size:20px;font-weight:700}.grid{stroke:#d9e2ec}.axis{stroke:#243b53;stroke-width:1.3}.label{fill:#52606d;font-size:12px}</style>",
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text class="title" x="70" y="34">{title}</text>',
    ]


def write_svg(path: Path, lines: list[str]) -> str:
    path.write_text("\n".join(lines + ["</svg>\n"]), encoding="utf-8")
    return str(path.relative_to(RESULT))


def x_pos(value: float, xmin: float, xmax: float, left: int, right: int) -> float:
    return left + (value - xmin) * (right - left) / (xmax - xmin)


def y_pos(value: float, ymin: float, ymax: float, top: int, bottom: int, log: bool) -> float:
    if log:
        value = math.log10(max(value, 1e-12))
        ymin = math.log10(max(ymin, 1e-12))
        ymax = math.log10(max(ymax, 1e-12))
    return bottom - (value - ymin) * (bottom - top) / (ymax - ymin)


def draw_axes(lines: list[str], ylabel: str, left: int, right: int, top: int, bottom: int) -> None:
    for i in range(5):
        y = top + i * (bottom - top) / 4
        lines.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}"/>')
    lines.append(f'<line class="axis" x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}"/>')
    lines.append(f'<line class="axis" x1="{left}" y1="{bottom}" x2="{left}" y2="{top}"/>')
    lines.append(f'<text class="label" transform="translate(24 310) rotate(-90)">{ylabel}</text>')


def draw_y_labels(
    lines: list[str],
    labels: list[tuple[str, float]],
    ymin: float,
    ymax: float,
    top: int,
    bottom: int,
    log: bool,
) -> None:
    for text, value in labels:
        y = y_pos(value, ymin, ymax, top, bottom, log)
        lines.append(f'<text class="label" x="42" y="{y + 4:.1f}">{text}</text>')


def draw_line_legend(
    lines: list[str],
    label: str,
    color: str,
    y: float,
    dashed: bool = False,
) -> None:
    dash = ' stroke-dasharray="6,5"' if dashed else ""
    lines.append(
        f'<line x1="690" y1="{y:.1f}" x2="724" y2="{y:.1f}" '
        f'stroke="{color}" stroke-width="3"{dash}/>'
    )
    lines.append(f'<circle cx="707" cy="{y:.1f}" r="3.5" fill="{color}"/>')
    lines.append(f'<text x="734" y="{y + 4:.1f}">{label}</text>')


def plot_iterations_vs_p(rows: list[dict[str, str]]) -> str:
    path = PLOTS / "iterations_vs_heterogeneity_ref5_np1.svg"
    left, right, top, bottom = 80, 650, 60, 430
    colors = ["#1d4ed8", "#dc2626", "#059669", "#7c3aed", "#ea580c"]
    p_values = sorted_values(rows, "p")
    ref = last_value(rows, "refinements")
    np = first_value(rows, "mpi_procs")
    selected = [
        row
        for row in rows
        if row["mpi_procs"] == np and row["refinements"] == ref
    ]
    max_iterations = positive_max([float(row["iterations"]) for row in selected], 10.0)
    y_max = 10 ** math.ceil(math.log10(max_iterations))
    lines = svg_start(f"iters vs p, ref={number_label(ref)}, np={number_label(np)}")
    draw_axes(lines, "CG iters, log", left, right, top, bottom)
    draw_y_labels(
        lines,
        [("1", 1), ("10", 10), ("100", 100), ("1k", 1000), ("10k", 10000)],
        1,
        y_max,
        top,
        bottom,
        True,
    )

    for index, prec in enumerate(PRECONDITIONERS):
        points = []
        for p_value in p_values:
            row = row_for(rows, np, p_value, ref, prec)
            if row:
                x = x_pos(float(p_value), float(p_values[0]), float(p_values[-1]), left, right)
                y = y_pos(float(row["iterations"]), 1, y_max, top, bottom, True)
                points.append((x, y))
        lines.append(
            f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in points)}" '
            f'fill="none" stroke="{colors[index]}" stroke-width="2.5"/>'
        )
        for x, y in points:
            lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{colors[index]}"/>')
        draw_line_legend(lines, prec, colors[index], 78 + 24 * index)

    for p_value in p_values:
        x = x_pos(float(p_value), float(p_values[0]), float(p_values[-1]), left, right)
        lines.append(f'<text class="label" x="{x-5:.1f}" y="454">{number_label(p_value)}</text>')
    lines.append('<text class="label" x="335" y="492">p</text>')
    return write_svg(path, lines)


def plot_total_time(rows: list[dict[str, str]]) -> str:
    path = PLOTS / "total_time_p5_ref5_np1.svg"
    left, right, top, bottom = 80, 760, 60, 430
    p = last_value(rows, "p")
    ref = last_value(rows, "refinements")
    np = first_value(rows, "mpi_procs")
    selected = [
        row
        for row in rows
        if row["mpi_procs"] == np and row["p"] == p and row["refinements"] == ref
    ]
    max_time = positive_max([float(row["total_s"]) for row in selected], 1.0)
    min_time = max(min(float(row["total_s"]) for row in selected), 1e-3)
    y_min = 10 ** math.floor(math.log10(min_time))
    y_max = 10 ** math.ceil(math.log10(max_time))
    lines = svg_start(
        f"total time, p={number_label(p)}, ref={number_label(ref)}, np={number_label(np)}"
    )
    draw_axes(lines, "seconds, log", left, right, top, bottom)
    draw_y_labels(
        lines,
        [(f"{value:g}", value) for value in [y_min, y_min * 10, y_min * 100, y_max] if value <= y_max],
        y_min,
        y_max,
        top,
        bottom,
        True,
    )
    slot = (right - left) / len(PRECONDITIONERS)
    colors = ["#1d4ed8", "#dc2626", "#059669", "#7c3aed", "#ea580c"]

    for index, prec in enumerate(PRECONDITIONERS):
        row = row_for(rows, np, p, ref, prec)
        if not row:
            continue
        value = float(row["total_s"])
        x = left + slot * index + slot * 0.22
        y = y_pos(value, y_min, y_max, top, bottom, True)
        lines.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{slot * 0.56:.1f}" height="{bottom-y:.1f}" fill="{colors[index]}"/>')
        lines.append(f'<text class="label" x="{x:.1f}" y="454">{prec}</text>')
        lines.append(f'<text class="label" x="{x:.1f}" y="{y-8:.1f}">{value:.3g}s</text>')
    return write_svg(path, lines)


def plot_strong_scaling(rows: list[dict[str, str]]) -> str:
    path = PLOTS / "strong_scaling_p5_ref5.svg"
    left, right, top, bottom = 80, 650, 60, 430
    colors = ["#1d4ed8", "#dc2626", "#059669", "#7c3aed", "#ea580c"]
    p = last_value(rows, "p")
    ref = last_value(rows, "refinements")
    np_values = sorted_values(rows, "mpi_procs")
    x_min = float(np_values[0])
    x_max = float(np_values[-1])
    speedups = []
    for prec in PRECONDITIONERS:
        serial = row_for(rows, np_values[0], p, ref, prec)
        if serial:
            serial_time = float(serial["total_s"])
            for np in np_values:
                row = row_for(rows, np, p, ref, prec)
                if row:
                    speedups.append(serial_time / float(row["total_s"]))
    y_max = max(x_max, math.ceil(positive_max(speedups, x_max)))
    lines = svg_start(f"strong scaling, p={number_label(p)}, ref={number_label(ref)}")
    draw_axes(lines, "speedup", left, right, top, bottom)
    draw_y_labels(
        lines,
        [(number_label(value), value) for value in range(0, int(y_max) + 1)],
        0,
        y_max,
        top,
        bottom,
        False,
    )
    ideal_points = [
        (x_pos(float(np), x_min, x_max, left, right), y_pos(float(np), 0, y_max, top, bottom, False))
        for np in np_values
    ]
    lines.append(
        f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in ideal_points)}" '
        f'fill="none" stroke="#64748b" stroke-width="2" stroke-dasharray="6,5"/>'
    )
    draw_line_legend(lines, "ideal", "#64748b", 198, dashed=True)

    for index, prec in enumerate(PRECONDITIONERS):
        serial = row_for(rows, np_values[0], p, ref, prec)
        if not serial:
            continue
        serial_time = float(serial["total_s"])
        points = []
        for np in np_values:
            row = row_for(rows, np, p, ref, prec)
            if row:
                speedup = serial_time / float(row["total_s"])
                points.append((x_pos(float(np), x_min, x_max, left, right), y_pos(speedup, 0, y_max, top, bottom, False)))
        lines.append(
            f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in points)}" '
            f'fill="none" stroke="{colors[index]}" stroke-width="2.5"/>'
        )
        for x, y in points:
            lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{colors[index]}"/>')
        draw_line_legend(lines, prec, colors[index], 78 + 24 * index)

    for np in np_values:
        x = x_pos(float(np), x_min, x_max, left, right)
        lines.append(f'<text class="label" x="{x-5:.1f}" y="454">{number_label(np)}</text>')
    lines.append('<text class="label" x="335" y="492">np</text>')
    return write_svg(path, lines)


def plot_setup_solve(rows: list[dict[str, str]]) -> str:
    path = PLOTS / "setup_solve_p5_ref5_np1.svg"
    left, right, top, bottom = 80, 720, 60, 430
    p = last_value(rows, "p")
    ref = last_value(rows, "refinements")
    np = first_value(rows, "mpi_procs")
    lines = svg_start(
        f"setup/solve, p={number_label(p)}, ref={number_label(ref)}, np={number_label(np)}"
    )
    draw_axes(lines, "seconds", left, right, top, bottom)
    methods = ["jacobi", "ssor", "ilu", "amg"]
    selected = [row_for(rows, np, p, ref, prec) for prec in methods]
    totals = [float(row["pc_setup_s"]) + float(row["solve_s"]) for row in selected if row]
    y_max = max(totals) * 1.2 if totals else 1.0
    y_labels = [(f"{y_max * fraction:g}", y_max * fraction) for fraction in [0, 0.25, 0.5, 0.75, 1]]
    draw_y_labels(lines, y_labels, 0, y_max, top, bottom, False)
    slot = (right - left) / len(methods)

    for index, prec in enumerate(methods):
        row = row_for(rows, np, p, ref, prec)
        if not row:
            continue
        setup = float(row["pc_setup_s"])
        solve = float(row["solve_s"])
        total = setup + solve
        setup_h = setup * (bottom - top) / y_max
        solve_h = solve * (bottom - top) / y_max
        x = left + slot * index + slot * 0.25
        width = slot * 0.50
        lines.append(f'<rect x="{x:.1f}" y="{bottom-solve_h:.1f}" width="{width:.1f}" height="{solve_h:.1f}" fill="#1d4ed8"/>')
        lines.append(f'<rect x="{x:.1f}" y="{bottom-solve_h-setup_h:.1f}" width="{width:.1f}" height="{setup_h:.1f}" fill="#ea580c"/>')
        lines.append(f'<text class="label" x="{x:.1f}" y="454">{prec}</text>')
        lines.append(f'<text class="label" x="{x:.1f}" y="{bottom-solve_h-setup_h-8:.1f}">{total:.3g}s</text>')

    lines.append('<rect x="750" y="72" width="14" height="14" fill="#ea580c"/><text x="772" y="84">setup</text>')
    lines.append('<rect x="750" y="98" width="14" height="14" fill="#1d4ed8"/><text x="772" y="110">solve</text>')
    return write_svg(path, lines)


def plot_optimality(rows: list[dict[str, str]]) -> str:
    """Iterations vs refinement level at fixed p, np=1. Flat lines = optimal."""
    path = PLOTS / "optimality_p5_np1.svg"
    left, right, top, bottom = 80, 650, 60, 430
    colors = ["#1d4ed8", "#dc2626", "#059669", "#7c3aed", "#ea580c"]
    p = last_value(rows, "p")
    np = first_value(rows, "mpi_procs")
    refs = sorted_values(rows, "refinements")
    if not refs:
        return ""

    selected = [
        row
        for row in rows
        if row["mpi_procs"] == np and row["p"] == p
    ]
    max_it = positive_max([float(row["iterations"]) for row in selected], 10.0)
    y_max = 10 ** math.ceil(math.log10(max_it))
    x_min, x_max = float(refs[0]), float(refs[-1])

    lines = svg_start(f"optimality, p={number_label(p)}, np={number_label(np)}")
    draw_axes(lines, "CG iters, log", left, right, top, bottom)
    draw_y_labels(
        lines,
        [(str(v), v) for v in [1, 10, 100, 1000, 10000] if v <= y_max],
        1, y_max, top, bottom, True,
    )

    for index, prec in enumerate(PRECONDITIONERS):
        points = []
        for ref in refs:
            row = row_for(rows, np, p, ref, prec)
            if row and float(row["iterations"]) > 0:
                x = x_pos(float(ref), x_min, x_max, left, right)
                y = y_pos(float(row["iterations"]), 1, y_max, top, bottom, True)
                points.append((x, y))
        if not points:
            continue
        lines.append(
            f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in points)}" '
            f'fill="none" stroke="{colors[index]}" stroke-width="2.5"/>'
        )
        for x, y in points:
            lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{colors[index]}"/>')
        draw_line_legend(lines, prec, colors[index], 78 + 24 * index)

    for ref in refs:
        x = x_pos(float(ref), x_min, x_max, left, right)
        lines.append(f'<text class="label" x="{x-5:.1f}" y="454">{number_label(ref)}</text>')
    lines.append('<text class="label" x="335" y="492">refinements</text>')
    return write_svg(path, lines)

    """Efficiency = t(base) / t(np) along (np, ref) pairs; ideal = 1.0."""
    path = PLOTS / "weak_scaling_p5.svg"
    left, right, top, bottom = 80, 650, 60, 430
    colors = ["#1d4ed8", "#dc2626", "#059669", "#7c3aed", "#ea580c"]
    p = last_value(rows, "p")
    if not p:
        return ""

    pairs = [(str(np), str(ref)) for np, ref in WEAK_SCALING_PAIRS]
    x_min, x_max = 1.0, float(WEAK_SCALING_PAIRS[-1][0])
    y_max = 1.2

    lines = svg_start(f"weak scaling, p={number_label(p)}")
    draw_axes(lines, "efficiency", left, right, top, bottom)
    draw_y_labels(
        lines,
        [(f"{v:g}", v) for v in [0, 0.25, 0.5, 0.75, 1.0]],
        0, y_max, top, bottom, False,
    )


    ideal = [
        (x_pos(float(np), x_min, x_max, left, right),
         y_pos(1.0, 0, y_max, top, bottom, False))
        for np, _ in pairs
    ]
    lines.append(
        f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in ideal)}" '
        f'fill="none" stroke="#64748b" stroke-width="2" stroke-dasharray="6,5"/>'
    )
    draw_line_legend(lines, "ideal", "#64748b", 198, dashed=True)

    for index, prec in enumerate(PRECONDITIONERS):
        base = row_for(rows, pairs[0][0], p, pairs[0][1], prec)
        if not base:
            continue
        base_time = float(base["total_s"])
        points = []
        for np, ref in pairs:
            row = row_for(rows, np, p, ref, prec)
            if row and float(row["total_s"]) > 0:
                eff = base_time / float(row["total_s"])
                points.append((x_pos(float(np), x_min, x_max, left, right),
                               y_pos(eff, 0, y_max, top, bottom, False)))
        if not points:
            continue
        lines.append(
            f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in points)}" '
            f'fill="none" stroke="{colors[index]}" stroke-width="2.5"/>'
        )
        for x, y in points:
            lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{colors[index]}"/>')
        draw_line_legend(lines, prec, colors[index], 78 + 24 * index)

    for np, _ in pairs:
        x = x_pos(float(np), x_min, x_max, left, right)
        lines.append(f'<text class="label" x="{x-5:.1f}" y="454">{np}</text>')
    lines.append('<text class="label" x="335" y="492">np</text>')
    return write_svg(path, lines)

def plot_iterations_vs_refinement(rows: list[dict[str, str]]) -> str:
    """Iterations vs mesh refinement for all p, np=1. Flat curves = h-optimal."""
    path = PLOTS / "iterations_vs_refinement_np1.svg"
    left, right, top, bottom = 80, 650, 60, 430
    colors = ["#1d4ed8", "#dc2626", "#059669", "#7c3aed", "#ea580c"]

    np = first_value(rows, "mpi_procs")
    refs = sorted_values(rows, "refinements")
    p_values = sorted_values(rows, "p")
    if not (refs and p_values):
        return ""

    p = p_values[-1]

    selected = [
        row for row in rows
        if row["mpi_procs"] == np and row["p"] == p
    ]
    max_it = positive_max([float(r["iterations"]) for r in selected], 10.0)
    y_max = 10 ** math.ceil(math.log10(max_it))
    x_min, x_max = float(refs[0]), float(refs[-1])

    lines = svg_start(
        f"iters vs refinement, p={number_label(p)}, np={number_label(np)}"
    )
    draw_axes(lines, "CG iters, log", left, right, top, bottom)
    draw_y_labels(
        lines,
        [(str(v), v) for v in [1, 10, 100, 1000, 10000] if v <= y_max],
        1, y_max, top, bottom, True,
    )

    for index, prec in enumerate(PRECONDITIONERS):
        points = []
        for ref in refs:
            row = row_for(rows, np, p, ref, prec)
            if row and float(row["iterations"]) > 0:
                x = x_pos(float(ref), x_min, x_max, left, right)
                y = y_pos(float(row["iterations"]), 1, y_max, top, bottom, True)
                points.append((x, y))
        if not points:
            continue
        lines.append(
            f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in points)}" '
            f'fill="none" stroke="{colors[index]}" stroke-width="2.5"/>'
        )
        for x, y in points:
            lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{colors[index]}"/>')
        draw_line_legend(lines, prec, colors[index], 78 + 24 * index)

    for ref in refs:
        x = x_pos(float(ref), x_min, x_max, left, right)
        lines.append(f'<text class="label" x="{x-5:.1f}" y="454">{number_label(ref)}</text>')
    lines.append('<text class="label" x="335" y="492">refinements</text>')
    return write_svg(path, lines)

def write_plots(rows: list[dict[str, str]]) -> list[str]:
    PLOTS.mkdir(parents=True, exist_ok=True)
    if not rows:
        return []
    return [
        plot_iterations_vs_p(rows),
        plot_total_time(rows),
        plot_strong_scaling(rows),
        plot_setup_solve(rows),
        plot_optimality(rows)
    ]


def main() -> int:
    RESULT.mkdir(parents=True, exist_ok=True)
    build_status = build_code()
    launcher = mpi_launcher()
    run_codes: dict[int, int] = {}
    rows: list[dict[str, str]] = []

    if build_status == 0:
        for np in MPI_PROCS:
            run_codes[np] = run_experiment(np, launcher)
            rows.extend(parse_run_log(np, run_codes[np]))

    write_csv(rows)
    plots = write_plots(rows)

    print(f"wrote {CSV_FILE}")
    for plot in plots:
        print(f"wrote {RESULT / plot}")

    return 0 if rows else 1


if __name__ == "__main__":
    raise SystemExit(main())
