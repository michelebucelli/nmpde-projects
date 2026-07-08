cmake_minimum_required(VERSION 3.13.4)
project(cardiac-electrophysiology LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(DEAL_II_DIR "/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/dealii/9.5.1/lib/cmake/deal.II/" CACHE PATH "Path to deal.II")

find_package(deal.II 9.0.0 QUIET CONFIG
  HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR})

if(NOT deal.II_FOUND)
  message(FATAL_ERROR
    "Could not find deal.II. Configure the project with -DDEAL_II_DIR=<path-to-deal.II> "
    "or add the installation prefix to CMAKE_PREFIX_PATH.")
endif()

message(STATUS "Using deal.II version ${DEAL_II_PACKAGE_VERSION} from ${DEAL_II_DIR}")

deal_ii_initialize_cached_variables()

if(NOT DEAL_II_WITH_MPI)
  message(FATAL_ERROR "Error! This project requires deal.II with MPI support.")
endif()

include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(cardiac_solver
  src/main.cpp
  src/CardiacProblem.cpp
)

deal_ii_setup_target(cardiac_solver)