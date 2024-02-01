# Preconditioning heterogeneous diffusion problems NMPDE Project
Implementation of preconditioning techniques for heterogeneous diffusion problems by Liborio Roman, Kamil Hanna and Filippo Mininno; for the Numerical Methods for Partial Differential Equations course by Prof.Quarteroni at Politecnico di Milano. Academic Year 2023-2024

## 1. Introduction
The objective of our project is to implement preconditioning techniques that aim to improve the computational efficiency (convergence rate) of solving heterogeneous diffusion problems, keeping in mind the optimality and parallel scalability of the problem.

## 2. Heterogeneous diffusion problems
Heterogeneous diffusion problems involve the spatially varying transport of a substance through a medium, where the diffusivity coefficient is not uniform across the domain. In mathematical terms, this is expressed by a partial differential equation of the form ∇ ⋅ (D(x) ∇u(x)) = f(x), where u(x) is the substance concentration, D(x) is the diffusivity coefficient, and f(x) represents any source or sink terms.

## 3. Our Implementation
Our project is structured as follows: 

The project is organized in the following folders:
- `src/`: Contains the files for both the serial & parallel implementations.
- `mesh/`: Contains a series of mesh files for testing.
- `results/`: Contains a csv file with the results of all the completed tests.
- `scripts/`: Contains a Python script for visualizing the results via plots.


## 4. Poisson problem
The 3D Poisson problem in which the diffusion coefficiant varies significantly (of orders of magnitude) over the domain:

$$ \begin{cases}
-\nabla \cdot (\mu \nabla u) = f \text{ in }\Omega \\
               u = 0 \text{ in } \partial \Omega \\
\end{cases} $$


Here, the coefficient function μ is defined as:

- μ(x) = 10^p if x is in region B
- μ(x)= 1 otherwise

$$
\[ \text{ where \( p > 0 \) } , B = \bigcup_{i=1}^{N} B_i, \text{ with } B_i \text{ a sphere of center } x_i \text{ and radius } r_i \]
$$

Weak Formulation :

If we multiply by a test function \(v\) and then integrate:
∫_{Ω} -∇ ⋅ (μ ∇u)v d𝑥 = ∫_{Ω} f v d𝑥 (1)

Considering the left-hand side:
∫_{Ω} -∇ ⋅ (μ ∇u)v d𝑥 = ∫_{Ω} μ ∇u ∇v d𝑥 - ∫_{∂Ω} μ v ∇u ⋅ n d𝑥 (2)

Where the last integral is zero due to boundary conditions.

Then:
∫_{Ω} μ ∇u ∇v d𝑥 = a(u, v)
= ∫_{Ω} f v d𝑥 = F(v) (3)


The problem in weak formulation now reads:
Find u in V: a(u, v) = F(v) for all v in V (4)

The problem in Galerkin formulation now reads:
Find u_h in V_h: a(u_h, v_h) = F(v_h) for all v_h in V_h (5)
where \(V_h = X_h^r \cap V\).


## 5. Dependencies and build
This project requires `CMake` for building the project,  `MPI` in order to enable the CPU parallelism.
In the case `openMP` is not installed or not found the project will still compile and run, only in fully serial manner.   

You can compile the project by running the following command in the root directory of the project:

```sh
mkdir build && cd build/ 
```
```sh
cmake ..
```
```sh
make
```
The following format

mpiexec -n (number of processors) ./Project_parallel -p (p-value) -P (preconditioner) -r (degree) -m (mesh file) -s (choose symmetry)

-p : set p value (10^p)
-P : set preconditionner
-r : set degree
-m : set mesh file (without path)
-s : Choose between symmetric (any value) /unsymmetric ("no")

Some examples...
```sh
./Project -p 4 -P amg -r 2 -m mesh-cube-40.msh -s 0
```

```sh
mpiexec -n 4 ./Project_parallel -p 4 -P amg -r 2 -m mesh-cube-40.msh -s 0
```
## 6. References
