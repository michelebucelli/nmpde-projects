cmake_minimum_required(VERSION 3.12.0)
project(matrix-free LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED "ON")

# Set default build type to Release.
if(NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
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

find_package(Kokkos 4.7.4 REQUIRED)

# Locate deal.II and initialize its variables.
find_package(deal.II 9.7.1 REQUIRED HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR})
deal_ii_initialize_cached_variables()

# Keep the Kokkos compile/include options, but remove only the static .a libraries from linkage
if(TARGET dealii::interface_kokkos)
  get_target_property(_libs dealii::interface_kokkos INTERFACE_LINK_LIBRARIES)
  list(FILTER _libs EXCLUDE REGEX "libkokkos")
  set_target_properties(dealii::interface_kokkos PROPERTIES INTERFACE_LINK_LIBRARIES "${_libs}")
endif()

# Add useful compiler flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLOG_TIMERS -Wall -Wfloat-conversion -Wmissing-braces -Wnon-virtual-dtor")

add_executable(main src/main.cpp)
deal_ii_setup_target(main)
