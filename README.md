# Preconditioning Heterogeneous Diffusion Problems

## Problem description

We solve the 3D Poisson problem with heterogeneous diffusion on Ω = (0, 1)³:

```
-∇·(μ ∇u) = f   in Ω
        u = 0   on ∂Ω
```

The diffusion coefficient μ(x) varies by orders of magnitude across the domain:
- μ(x) = 10^p if x lies inside any of 12 predefined spheres
- μ(x) = 1 otherwise
  where p > 0 controls the strength of heterogeneity.

We investigate the effectiveness, optimality, and parallel scalability of five
preconditioning strategies available in the deal.II `TrilinosWrappers` namespace

The provided `setup_modules.sh` script handles this for you. It must be `source`d for the changes to take effect in your current shell.

## Requirements

- deal.II ≥ 9.5.1 (with Trilinos and MPI support)
- CMake ≥ 3.13
- C++17 compiler
- Python 3.8+ (for the benchmark script; no external packages needed)


```bash
source setup_modules.sh
```

## Building

```bash
mkdir -p build && cd build
cmake ..
make -j$(nproc)
cd ..
```

### 1. Verification mode

Runs a manufactured-solution convergence test with μ = 1 (homogeneous),
using u(x,y,z) = sin(πx) sin(πy) sin(πz). Prints an L²/H¹ convergence
table. Expected rates: L² = 2.0, H¹ = 1.0 (continuous P1 elements on
tetrahedra).

```bash
./build/elliptic verify              # refinements 3 to 5
./build/elliptic verify 2 6          # refinements 2 to 6
mpirun -n 4 ./build/elliptic verify  # parallel
```

### 2. Single-case mode

Run one (p, refinement, preconditioner) combination. Writes `.pvtu` output
for ParaView.

```bash
# All 5 preconditioners at p=3, refinement 4
./build/elliptic 3 4
 
# Only AMG at p=5, refinement 5, with 4 MPI processes
mpirun -n 4 ./build/elliptic 5 5 amg
```

### 3. Full sweep mode (no arguments)

Sweeps all combinations of p ∈ {1,2,3,4,5}, refinements ∈ {3,4,5}, and
all 5 preconditioners. Used by the benchmark script.

```bash
./build/elliptic
```

## Running the full benchmark

The Python script automates the sweep across MPI process counts {1, 2, 3, 4},
collects structured logs, and produces:

```bash
python3 src/script/benchmark_preconditioners.py
```

**Outputs:**
- `result/benchmark_results.csv` — All timing and iteration data
- `result/plots/iterations_vs_heterogeneity_ref5_np1.svg` — CG iterations vs p
- `result/plots/total_time_p5_ref5_np1.svg` — Total solve time comparison
- `result/plots/strong_scaling_p5_ref5.svg` — Speedup vs MPI processes
- `result/plots/setup_solve_p5_ref5_np1.svg` — Setup vs solve time breakdown
- `result/plots/optimality_p5_np1.svg` — Iterations vs mesh refinement (h-optimality)
## Visualization

ParaView output (`.pvtu` files) is written to `output/` and includes:
- `solution` — the computed solution field
- `mu` — the diffusion coefficient field (showing the sphere inclusions)
```bash
# Generate output for a specific case
mpirun -n 4 ./build/elliptic 5 4 amg
# Open output/heterogeneous-amg-p5-ref4.pvtu in ParaView
```

## Report and presentation

The report in `report/report.tex` explains the strong, weak, and discrete
formulations, the CG/preconditioner algorithms, the distributed-memory
implementation, and the benchmark results. Rebuild it with `make report`.
The short presentation outline is in `presentation.md` and follows the same
notation and results.

## Verification results

```
--- convergence verification (mu = 1) ---
    h             L2               H1
1.2500e-01 5.7462e-03     - 2.1818e-01     -
6.2500e-02 1.4367e-03  2.00 1.0905e-01  1.00
3.1250e-02 3.5919e-04  2.00 5.4524e-02  1.00
```
