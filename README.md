# WaveFEM Solver

deal.II solver for Project 2, Numerical Methods for Partial Differential
Equations, Politecnico di Milano. The repository contains a sequential baseline
and a separate MPI implementation.

## Team

- Salvatore Mariano Librici - `salvatoremariano.librici@mail.polimi.it`
- Rong Huang - `rong.huang@mail.polimi.it`
- Hirdesh Kumar - `hirdesh.kumar@mail.polimi.it`
- Mehdi Ghiasipour - `mehdi.ghiasipour@mail.polimi.it`

## What the Code Does

It solves $u_{tt}-\Delta u=f$ on $(0,1)^2$, with $u=0$ on the boundary,
initial displacement $\sin(\pi x)\sin(\pi y)$, and zero initial velocity.

The mesh contains **squares, not triangles**: `subdivided_hyper_cube`
creates $N^2$ cells and `FE_Q(1)` gives continuous bilinear Q1 elements.
The unknowns are nodal values. Each square contributes local mass/stiffness
matrices and a load vector, assembled using its global DoF indices.
Degree p selects Qp; the reported convergence tests use p=1.

The weak form is $(u_{tt},v)+(\nabla u,\nabla v)=(f,v)$ for $v\in H^1_0$.
It gives $M U''+K U=F$. The mesh and unit coefficients are fixed, so M and K
are assembled once. Newmark average acceleration uses beta=1/4, gamma=1/2:
each step solves $(M+\Delta t^2 K/4)A^{n+1}=F^{n+1}-K U_{\rm pred}$, then
updates displacement and velocity. This matrix is also constant. No direct
factorization or cached preconditioner is used.

CG+SSOR is the default. CG is valid because the consistent mass matrix is
positive definite, stiffness is nonnegative, and boundary **rows and columns**
are eliminated symmetrically. GMRES and Jacobi are configurable alternatives.
See the [measured comparison](submission/assets/tables/solver-comparison.md);
the single sine mode is a special benchmark, not evidence of general optimality.

## Build and Run

Inside the course environment with deal.II 9.5.1:

```bash
bash scripts/docker_build_run.sh
bash scripts/docker_build_run.sh 32 1 1.0 0.0025 10 output_forced forced cg jacobi
```

Defaults: N=24, degree=1, T=1, dt=0.005, output every 10 steps, free, CG+SSOR.
T/dt must be a positive integer. The helper loads course modules when available
and reuses the build. On another deal.II installation, configure its CMake
search path yourself; the course module names are not a Docker image installer.

## Validation and Results

The exact solution is $u=\sin(\pi x)\sin(\pi y)\cos(\omega t)$, giving
$f=(2\pi^2-\omega^2)u$. Cases: free uses $\omega=\sqrt{2}\pi$;
forced uses $\omega=\pi$.

```bash
source scripts/common.sh
load_course_modules
python scripts/review_validation.py
```

This runs spatial refinement and halved-dt sensitivity checks, all four
solver/preconditioner combinations, sparse matrix checks, a Q2 smoke case,
invalid-input tests, and VTU field/boundary checks. It writes reproducible
[convergence](submission/assets/tables/convergence.md),
[energy](submission/assets/tables/energy-summary.md), and comparison tables.
Raw outputs/logs go to ignored `output/validation/review/`.

Errors are **absolute displacement** L2 norms and H1 seminorms integrated with
3x3 Gauss points for Q1. Error evolution versus time uses semilog-y axes.
Separate convergence plots show final errors at T=1 versus
h=sqrt(2)/N on log-log axes. Smooth Q1 approximation predicts orders 2 and 1;
the experiment checks temporal contamination before interpreting those rates.

Energy is $E_h=(V^TMV+U^TKU)/2$. Free Newmark energy is conserved in exact
arithmetic; the table measures numerical drift. Forced energy varies because
the load does work. A conserved energy does not prove phase accuracy.

The optional `probe.csv` records the displacement at `(x,y)=(0.5,0.5)` and is
used by the dispersion check to compare numerical and exact peak times.

## What to Open

- `output/wave_solution.pvd` in ParaView: select **u**, then **Surface With Edges**
  or **Warp By Scalar**. Height represents displacement; the PDE remains 2D.
- `output/diagnostics.csv`, `errors.png`, `energy.png`, `energy_change.png`.
- [Report PDF](submission/final-report/latex/WaveFEM_Solver_Report.pdf).
- [Oral presentation notes](PRESENTATION_NOTES.md).
- [Technical review and remaining concerns](TECHNICAL_REVIEW.md).
- [Commands](COMMANDS.md), including every figure/table reproduction command.

## What to Study

Use this order for the oral presentation:

1. `PRESENTATION_NOTES.md`: short answers to likely questions.
2. This README: the implemented PDE, discretization and time-stepping chain.
3. `src/sequential/WaveEquation.cpp`: baseline setup, assembly and Newmark loop.
4. `src/mpi/WaveEquationMPI.cpp`: the same problem with distributed mesh and vectors.
5. `submission/assets/tables/`: measured convergence, energy, solver and MPI results.
6. `submission/final-report/latex/WaveFEM_Solver_Report.pdf`: complete report.

The final 18-slide deck is available at
`presentation/WaveFEM_Solver_Presentation.pptx`.

The main chain is:
`PDE -> weak form -> Q1 FEM -> M/K -> Newmark -> SPD solve -> verification -> output`.

## HPC Scope

The default `wave_fem` target is serial. The separate `wave_fem_mpi` target is
a completed lightweight distributed-memory implementation using a distributed
mesh, Trilinos vectors/matrices, MPI reductions and distributed VTK output.
The serial target remains the reference for controlled numerical verification;
the MPI target is the HPC execution path.

Run the MPI target on a machine where multi-rank Open MPI is configured:

```bash
source scripts/common.sh
load_course_modules
cmake --build build -j2
mpirun -np 4 ./build/wave_fem_mpi 64 0.05 0.005 output/mpi-scaling/manual-np4 free
```

The MPI run writes a final `.pvtu` snapshot per run. The scaling script stores
results in `output/mpi-scaling/`. `run_scalability_study.sh`
is still a serial problem-size study; it is not strong or weak MPI scaling.
The reported comparison is part of the HPC study carried out for this project.
For the same `N=192`, `T=0.1`, `dt=0.0025` problem, the serial executable took
14.533 s, while MPI took 13.823 s with 1 rank, 11.529 s with 2 ranks and
7.339 s with 4 ranks. Final numerical errors agree. MPI-4 gives speedup 1.98
and 49.5% parallel efficiency relative to serial. The MPI time loop decreases
from 12.322 s to 6.361 s; the end-to-end comparison uses different sparse
backends in the two targets.

Source is in `src/`, scripts in `scripts/`, and report material in
`submission/`. Build and raw solver outputs are not committed.
