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


# CMake shenanigans
From a fresh CMake config, the project doesn't compile because it misses the `bench_config.h` file. From previous discussions that file is generated at config time by cmake based on config variables, but there is no mention of the generation in the root `CMakeLists.txt`.

This is introducing some strange coupling between the benchmark infraastructure and the matrix free library. I propose to eliminate this complexity


# SIMD_WIDTH
We want the simd width to be specified by the user (why?)
In any case, if we don't specify it deal.ii will provide the maximum width supported by the target.

The case we want to absolutely avoid is when the user selects at compile time a config
