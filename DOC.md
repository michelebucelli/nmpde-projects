# Parametrization
The `MatrixFreeADR` class is parametric on:

| name          | type        | description                                     |
|---------------|-------------|-------------------------------------------------|
| dim           | int         | the dimension of the problem                    |
| fe_degree     | int         | the polynomial degree for the FE discretization |
| number        | categorical | which floating point type to use                |
| simd_width    | int         | how big should the vector register be           |
| AdvectionFunc | func        | advection coefficient                           |
| DiffusionFunc | func        | diffusion coefficient                           |
| ReactionFunc  | func        | reaction coefficient                            |

we encode those parameters as template parameters of the `MatrixFreeADR` class. This choice is obligated for `dim`, `fe_degree`, `number` and `simd_width` as we are using deal.ii templated classes.

We use the functor pattern for passing problem coefficients to the `MatrixFreeADR` class. If we were to pass function objects, the compiler would have a hard time inlining the coefficient functions when computing the matrix action in `MatrixFreeADR<>::solve(...)`. A matrix free solver evaluates the coefficients every time as opposed to a matrix based solver, so the cost of inefficiencies could quickly compound.

By passing the functions as template parameters, the `MatrixFreeADR` statically owns the implementation of those functions, this helps the compiler inline the calls. We have not measured what is the performance benefit of our approach.

## Discussion on parametrizations
  * the claim of templating on the coefficient functions must be proved with an experiment, but it's an ass experiment
  * do we really want to template on `simd_width`? If we don't specify the dimension, deal.ii will automatically give the best width. Maybe we can measure the throughput when varying the simd width, but it makes sense that the best solution is having the widest register supported by the target architecture, am I missing something?
