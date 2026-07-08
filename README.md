## matrix-free operators in deal.II

- Luca Ballone 10813225 (lechebs)
- Matteo Mei 10706418 (FattiMei)
- Tommaso Angelaccio 11089177 (tiker02)
- Fabio Ghattas 10764040 (FabioGhattasPolimi)

### Compiling the program

The project was developed and tested on deal.II 9.5.1. and 9.7.1.

Compile using CMake:
```
mkdir build && cd build
cmake ..
make
```

Note that the compilation of the targets `single-node` and `tube_poiseuille` may take some time, since the build compiles the matrix-free implementation for many `fe_degree` values. To build and execute only the convergence test, run:
```
make main
./main
```

JIT compilation to support runtime `fe_degree` values, as well as the patch to deal.II headers to improve GPU performance under Kokkos, can be found at https://github.com/lechebs/mf-adr/tree/jit.

## Portable matrix-free operators in deal.II

### Installing deal.II with Kokkos for CUDA support (Ubuntu)

Install the required dependencies with:
```
sudo apt install libopenmpi-dev libp4est-dev
```

Build Kokkos 4.7.4 from source with CUDA and OpenMP backends:
```
wget https://github.com/kokkos/kokkos/releases/download/4.7.04/kokkos-4.7.04.tar.gz
tar -xvf kokkos-4.7.04.tar.gz && cd kokkos-4.7.04
cmake -B build -DCMAKE_CXX_STANDARD=20 -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_<MICROARCHITECTURE><COMPUTE_CAPABILITY>=ON -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_DEPRECATION_WARNINGS=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON -DKokkos_ENABLE_CUDA_CONSTEXPR=ON -DKokkos_ENABLE_DEPRECATED_CODE_4=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON
cmake --build build -jN
mkdir <PATH_TO_KOKKOS_INSTALL_DIR>
cmake --install build --prefix <PATH_TO_KOKKOS_INSTALL_DIR>
```

More information on configuring the Kokkos build can be found at https://kokkos.org/kokkos-core-wiki/get-started/building-from-source.html.

Build deal.II 9.7.1 from source with MPI, Kokkos and p4est support with:
```
wget https://github.com/dealii/dealii/releases/download/v9.7.1/dealii-9.7.1.tar.gz
tar -xvf dealii-9.7.1.tar.gz && cd dealii-9.7.1
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=<PATH_TO_DEALII_INSTALL_DIR> -DCMAKE_CXX_COMPILER=<PATH_TO_KOKKOS_INSTALL_DIR>/bin/nvcc_wrapper -DCMAKE_CXX_STANDARD=20 -DDEAL_II_WITH_KOKKOS=ON -DDEAL_II_WITH_MPI=ON -DDEAL_II_WITH_P4EST=ON -DDEAL_II_MPI_WITH_DEVICE_SUPPORT=ON -DKOKKOS_DIR=<PATH_TO_KOKKOS_INSTALL_DIR> ..
make -jN install
```

More details on configuring the deal.II build can be found at https://dealii.org/developer/readme.html.

### Compiling the program

The program can be compiled with:
```
mkdir build && cd build
export CMAKE_PREFIX_PATH=/path/to/dealii/install/dir
cmake -DBUILD_WITH_KOKKOS=ON ..
make
```

The compilation targets and executable produced are the same as the CPU version.

### References

- https://dealii.org/current/doxygen/deal.II/step_37.html
- https://dealii.org/current/doxygen/deal.II/step_64.html
- https://dealii.org/current/doxygen/deal.II/group__matrixfree.html
- https://dealii.org/current/doxygen/deal.II/classMatrixFreeOperators_1_1Base.html
- https://dealii.org/current/doxygen/deal.II/classMatrixFree.html
- https://dealii.org/current/doxygen/deal.II/classFEEvaluation.html


