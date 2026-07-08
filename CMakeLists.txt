cmake_minimum_required(VERSION 3.12.0)
project(matrix-free LANGUAGES CXX C)

option(BUILD_WITH_KOKKOS "Build the Kokkos-based portable implementation" OFF)

set(CMAKE_CXX_STANDARD 17)
if (BUILD_WITH_KOKKOS)
  set(CMAKE_CXX_STANDARD 20)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED "ON")

# Set default build type to Release.
if(NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()
message(STATUS)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  add_definitions(-DBUILD_TYPE_DEBUG)
endif()

# Locate MPI compiler.
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")

set(DEAL_II_MIN_VERSION 9.3.1)

if (BUILD_WITH_KOKKOS)
  find_package(Kokkos 4.7.4)
  set(DEAL_II_MIN_VERSION 9.7.1)
endif()

# Locate deal.II and initialize its variables.
find_package(deal.II ${DEAL_II_MIN_VERSION} REQUIRED HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR})
deal_ii_initialize_cached_variables()

# search for likwid, but make it quiet because
# not all users have installed it
find_package(likwid QUIET)
if (likwid_FOUND)
  add_definitions(-DUSE_LIKWID -DLIKWID_PERFMON)
  message(STATUS "Found likwid")
endif()

# Add useful compiler flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfloat-conversion -Wmissing-braces -Wnon-virtual-dtor")

if (BUILD_WITH_KOKKOS)
  add_compile_definitions(BUILD_WITH_KOKKOS)

  # To be used if kokkos was linked statically against deal.II
  if(TARGET dealii::interface_kokkos)
    get_target_property(_libs dealii::interface_kokkos INTERFACE_LINK_LIBRARIES)
    list(FILTER _libs EXCLUDE REGEX "libkokkos")
    set_target_properties(dealii::interface_kokkos PROPERTIES INTERFACE_LINK_LIBRARIES "${_libs}")
  endif()
endif()

add_executable(main src/main.cpp)
target_link_libraries(main)
target_compile_definitions(main PUBLIC)
deal_ii_setup_target(main)

add_executable(single-node benchmark/single_node.cpp)
target_include_directories(single-node PRIVATE "src/")
target_link_libraries(single-node)
deal_ii_setup_target(single-node)

if (likwid_FOUND)
  target_link_libraries(single-node likwid::likwid)
endif()

add_executable(tube_poiseuille tube_poiseuille/experiment.cpp)
target_include_directories(tube_poiseuille PRIVATE "src/")
target_link_libraries(tube_poiseuille)
deal_ii_setup_target(tube_poiseuille)
