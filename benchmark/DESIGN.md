# deal.ii versioning
Ubuntu distribution packages ship with version 9.5.1. Version 9.7.1 was required in the `kokkos` branch, but it's incompatible when defining affine constraints inside the matrix free code.

Here is the error from compiling the benchmark application
```
/home/matteo/github/mf-adr/src/MatrixFreeADR.hpp|118 col 23| error: no matching function for call to ‘dealii::AffineConstraints<double>::reinit(const dealii::IndexSet&, const dealii::IndexSet&)’
```

After a team discussion we decided to revert the changes such that the project can be compiled with version 9.5.1 and possibly add a compile time selection between the two implementations. The modifications would be verified by a passing convergence test on the main target


```cpp
        const IndexSet locally_relevant_dofs =
            DoFTools::extract_locally_relevant_dofs(dof_handler);

        constraints.reinit(locally_relevant_dofs);
```
