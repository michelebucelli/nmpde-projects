cmake_minimum_required(VERSION 3.12.0)

# Use toolchain's MPI wrappers to avoid glibc conflicts
set(CMAKE_C_COMPILER "/u/sw/toolchains/gcc-glibc/11.2.0/base/bin/mpicc")
set(CMAKE_CXX_COMPILER "/u/sw/toolchains/gcc-glibc/11.2.0/base/bin/mpicxx")

project(elliptic LANGUAGES CXX C)

# Export a compilation database for clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set the path to the Boost library
set(Boost_DIR "/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/boost/1.76.0/lib/cmake/Boost-1.76.0")

# Set the path to the deal.II library
set(DEAL_II_DIR "/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/dealii/9.5.1/lib/cmake/deal.II")

include(common/cmake-common.cmake)

add_executable(elliptic
  src/src_elliptic/main.cpp
  src/src_elliptic/HeterogeneousDiffusion.cpp)
  
deal_ii_setup_target(elliptic)
