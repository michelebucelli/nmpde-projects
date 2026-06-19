How to benchmark a dealii application

## Generalities
On a single node I want to measure the performance of a given iterative method.
I'm not concerned with the convergence properties of the method, just on its implementation.
No preconditioner is used (`PreconditionerIdentity`) when comparing matrix based and matrix free implementations.

I'm proposing to measure the *time per iteration / dofs* by running a fixed amount of iterations of the method while sweeping the degrees of freedom (independent variable).
Taking the average of a 100 runs is a technique to reduce the variance of the measurement since fast iterations are more exposed to fluctuations.

## Single node scaling
There are two methods for controlling the dofs:
  * **h-scaling:** increases the resolution of the mesh
  * **p-scaling:** increases the number of dofs per cell

both are worth exploring since they expose different aspects of the computing architecture.

### p-scaling
If there are many more dofs per cell it means that the cost of computing the geometric transformations from the reference cell to the real cell is amortized.
This could improve the overall iteration time and the arithmetic intensity.
The effect should be more significant for unstructured meshes.

**TODO:** p-scaling is problematic when the solver code is templated on `fe_degree`. How can we solve this? In principle we don't want to recompile the application, maybe we can do some dispatching directly from the code.

### Unstructured meshes
Challenges arise from the use of unstructured meshes:
  * each cell may require its own geometrical transformation
  * the connectivity of the mesh is irregular
  * external programs (`gmsh`) must be used to generate such meshes


## Experiment 1: single node hp-scaling on structured meshes
```
for every method:
    for many mesh sizes:
        for many degrees:
            profile 100 runs of the iterative method
```

We will measure:
  * time per iteration
  * arithmetic intensity
  * number of vectorized instructions
  * memory bandwidth

Produce an executable which is parametric on the triplet (method, mesh size, fe_degree) and use external scripts to collect the results.
Python is preferred instead of bash because of the plotting ecosytem and the need of implementing sweeping logic.
Refer to the `subprocess` library.


## Experiment 2: single node p-scaling on structured vs. unstructured mesh
```
for every method:
    for every mesh type:
        for many degrees:
            profile 100 runs of the iterative method
```

Given the difficulty of generation, we could use an unstructured mesh from Bucelli and compare it with a structured mesh of the same cell number.
If an equivalent structured mesh is not present, we can generate one with `gmsh`.
We don't scale with the number of cells in this case, this simplifies the mesh generation.

Note that we may still use the same executable of experiment 1. We just need to come up with a proper command line interface.
Maybe the mesh type and mesh size could be bound to the same command line argument. If the argument is an integer, we know we must generate the mesh, otherwise we load it from file.

