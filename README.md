## Portable matrix-free operators in deal.II

### Installing deal.II with Kokkos (Ubuntu)

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
cmake ..
make
```

### Useful resources

- https://dealii.org/current/doxygen/deal.II/step_37.html
- https://dealii.org/current/doxygen/deal.II/step_64.html

- https://dealii.org/current/doxygen/deal.II/group__matrixfree.html
- https://dealii.org/current/doxygen/deal.II/classMatrixFreeOperators_1_1Base.html
- https://dealii.org/current/doxygen/deal.II/classMatrixFree.html
- https://dealii.org/current/doxygen/deal.II/classFEEvaluation.html

### TODO

- [ ] various TODOs in the code
- [ ] generic preconditioners
- [ ] multigrid
- [ ] benchmark pipeline
- [ ] convergence test suite
=======
