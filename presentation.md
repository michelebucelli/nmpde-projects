# Preconditioning Heterogeneous Diffusion Problems

## Slide 1 — Problem and objective

- Solve `-div(mu grad u) = 1` in `(0,1)^3`, with `u=0` on the boundary.
- Twelve spherical inclusions have `mu=10^p`; the background has `mu=1`.
- Compare identity, Jacobi, SSOR, ILU, and AMG with distributed CG.

## Slide 2 — Weak and discrete problems

- `V = H_0^1(Omega)` and `a(u,v) = integral mu grad(u) dot grad(v)`.
- Find `u` such that `a(u,v) = integral f v` for every `v` in `V`.
- Use continuous `P1` finite elements on tetrahedra:
  `A_ij = integral mu grad(phi_j) dot grad(phi_i)`, `b_i = integral f phi_i`.
- The unfitted spherical interfaces are sampled at quadrature points.

## Slide 3 — Solver and implementation

- CG starts from zero and stops at `||r|| <= 10^-10 ||b||`.
- Distributed deal.II triangulation, DoFs, sparse matrix, vectors, and MPI reductions.
- AMG supplies multilevel coarse correction; Jacobi/SSOR/ILU provide increasingly local algebraic
  information.
- Record preconditioner setup, solve time, iterations, and CG/Lanczos condition estimate.

## Slide 4 — Verification

- Manufactured solution: `u = sin(pi x) sin(pi y) sin(pi z)`, `mu=1`.
- Observed under mesh halving: `L2` rate `2.00`, `H1` rate `1.00`.
- This validates assembly, constraints, distributed solve, and error integration independently of
  heterogeneous-interface accuracy.

## Slide 5 — Robustness to contrast

- Finest serial mesh, `p=5`: identity needs 16,067 iterations; AMG needs 31.
- Estimated condition number falls from `1.05e7` (identity) to `18.5` (AMG).
- Local methods help, but their iteration counts still grow strongly with contrast.

## Slide 6 — Cost and parallelism

- Hardest serial total time: identity `16.941 s`, AMG `0.262 s`.
- At four MPI ranks AMG remains fastest at `0.168 s`; speedup is `1.57x`.
- Small systems make launch, communication, and coarse-grid costs visible; speedup and absolute
  time answer different questions.

## Slide 7 — Conclusions

- AMG is the best default in the tested heterogeneous elliptic regime.
- Jacobi is attractive for cheap, mild cases; ILU is a useful intermediate option.
- Stronger claims require repeated timings and larger meshes; the repository contains the source,
  benchmark driver, raw CSV, plots, and reproducible report.
