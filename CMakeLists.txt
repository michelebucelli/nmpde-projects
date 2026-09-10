cmake_minimum_required(VERSION 3.12.0)
project(WaveFEMSolver LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

message(STATUS)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS)

# The course laboratories build deal.II programs through MPI compilers, even for
# serial examples. We keep the same convention so the project behaves like the
# lab code on the official Docker/HPC environment.
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")

find_package(Boost 1.72.0 REQUIRED
  COMPONENTS filesystem iostreams serialization
  HINTS ${BOOST_DIR} $ENV{BOOST_DIR} $ENV{mkBoostPrefix})
include_directories(${Boost_INCLUDE_DIRS})

find_package(deal.II 9.3.1 REQUIRED
  HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR} $ENV{mkDealiiPrefix})
deal_ii_initialize_cached_variables()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfloat-conversion -Wmissing-braces -Wnon-virtual-dtor")

add_executable(wave_fem
  src/sequential/main.cpp
  src/sequential/WaveEquation.cpp
)
deal_ii_setup_target(wave_fem)

add_executable(wave_fem_mpi
  src/mpi/main_mpi.cpp
  src/mpi/WaveEquationMPI.cpp
)
deal_ii_setup_target(wave_fem_mpi)
