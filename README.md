# Hybrid Matrix-Free Solver

This project implements a hybrid MPI + multithreading matrix-free finite element
solver for the convection–diffusion–reaction problem

$$
\begin{cases}
-\nabla \cdot (\mu \nabla u) + \boldsymbol{\beta} \cdot \nabla u + \gamma u = f,
  & \text{in } \Omega, \\
u = g, & \text{on } \Gamma_D \subset \partial\Omega, \\
\nabla u \cdot \mathbf{n} = h, & \text{on } \Gamma_N = \partial\Omega \setminus \Gamma_D.
\end{cases}
$$

## Building

```bash
$ mkdir -p build
$ cd build
$ cmake ..
$ make
```

## Running

The solver is configured by a deal.II parameter file. Write the fully
documented default one, then edit it:

```bash
$ ./build/hybrid_matrix_free_solver --generate input.prm
$ ./build/hybrid_matrix_free_solver -i input.prm
```

The file covers the problem (`Dimension`, `Degree`, `Mu`, `Gamma`, `Beta`),
the discretization (`Refinements`, `Cycles`), the solver (`Preconditioner`,
`Tolerance`), the TBB thread count and the output. The MPI rank count is not
a parameter; it comes from `mpirun -n`.

Single entries can be overridden on the command line, so a sweep varies one
knob without a file per data point:

```bash
$ mpirun -n 4 ./build/hybrid_matrix_free_solver -i input.prm \
      --threads 2 --set Solver/Preconditioner=jacobi
```

`--threads N` and `--refine N` are shorthands for `--set
Parallelism/Threads=N` and `--set Discretization/Refinements=N`. Arguments the
solver does not recognise are forwarded to MPI. `--help` lists everything.

Every run prints a machine readable `SUMMARY` line with the rank and thread
counts, the DoF count, the iteration count and the setup, rhs and solve wall
times. Setting `Output/Benchmark applications` to a non-zero value adds a
`VMULT` line timing bare operator applications, isolated from the Krylov and
multigrid work.

## Benchmarks

`scripts/benchmark.sh` runs the four suites of the study and writes one CSV per
suite into `results/`: convergence of both backends, matrix-based against
matrix-free, MPI strong scaling, and the hybrid rank/thread grid.

`Solver/Backend` picks the operator representation: `matrix_free` or
`matrix_based`, the latter assembling a Trilinos sparse matrix from the same
bilinear form as the reference for the comparison.

## Testing

```bash
$ cd build && ctest
```

## Repository layout

```
include/    cdr_problem     problem data and manufactured solution
            cdr_operator    matrix-free CDR operator
            mg_preconditioner  geometric multigrid V-cycle
            cdr_solver      the FGMRES driver and the convergence table
            cdr_parameters  the parameter file
            cdr_assembled   matrix-based reference solver
            cli             the command line
src/        main.cpp
tests/
scripts/    benchmark.sh    the benchmark sweeps
results/    *.csv           measured data
doc/        report.pdf      the report
input.prm   generated default parameter file
```
