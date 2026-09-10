# MPI scaling measurement

Same Q1 free-wave problem for every run: `N=192`, `T=0.1`, `dt=0.0025`.
The table compares the original serial executable with the MPI executable.
MPI phase values are the maximum across ranks; serial has no MPI phase split.

| Implementation | Ranks | DoFs | Cells | Setup (s) | Time loop (s) | Post-process (s) | Total (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| serial | 1 | 37249 | 36864 | -- | -- | -- | 14.533 |
| MPI | 1 | 37249 | 36864 | 1.177 | 12.322 | 0.324 | 13.823 |
| MPI | 2 | 37249 | 36864 | 1.148 | 10.183 | 0.198 | 11.529 |
| MPI | 4 | 37249 | 36864 | 0.865 | 6.361 | 0.158 | 7.339 |

All implementations run correctly and produce matching final errors. With the
same diagnostics and output policy, MPI-1 is slightly faster than the serial
baseline, and MPI-4 reaches a speedup of 1.98 with an efficiency of 49.5%
relative to serial. The MPI time-loop breakdown decreases from 12.322 s on one
rank to 6.361 s on four ranks. The efficiency is below ideal because each rank
still participates in synchronisation and distributed sparse solves.
