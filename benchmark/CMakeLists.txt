cmake_minimum_required(VERSION 3.12.0)
project(matrix-free LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
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

# Locate Boost.
find_package(Boost 1.72.0 REQUIRED
  COMPONENTS filesystem iostreams serialization
  HINTS ${BOOST_DIR} $ENV{BOOST_DIR} $ENV{mkBoostPrefix})
message(STATUS "Using the Boost-${Boost_VERSION} configuration found at ${Boost_DIR}")
message(STATUS)
include_directories(${Boost_INCLUDE_DIRS})

# Locate deal.II and initialize its variables.
find_package(deal.II 9.3.1 REQUIRED
  HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR} $ENV{mkDealiiPrefix})
deal_ii_initialize_cached_variables()

# Add useful compiler flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLOG_TIMERS -Wfloat-conversion -Wmissing-braces -Wnon-virtual-dtor -mfma")

#benchmark compiler flags
set(FE_DEGREE 1 CACHE STRING "Polynomial degree of the Finite Elements used")
set(PROBLEM_DIMENSIONS 2 CACHE STRING "Mesh space dimensions")
set(SIMD_WIDTH 4 CACHE STRING "Size param for deal.II VectorizedArray")
set(NUMBER_TYPE "double" CACHE STRING "Floating point type (float, double)")

configure_file(
"${CMAKE_CURRENT_SOURCE_DIR}/benchmark_config.h.in"
"${CMAKE_CURRENT_BINARY_DIR}/bench_config.h"
)

file(GLOB SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../src/*.hpp")
add_library(mf-adr INTERFACE)
target_include_directories(mf-adr INTERFACE  "${CMAKE_CURRENT_SOURCE_DIR}/../src" "${CMAKE_CURRENT_BINARY_DIR}")


set(TARGET_NAME "bex_d${PROBLEM_DIMENSIONS}_p${FE_DEGREE}_w${SIMD_WIDTH}_${NUMBER_TYPE}")


message(STATUS "Building mf simulation in ${PROBLEM_DIMENSIONS} dimensions with poynomial degree ${FE_DEGREE}, using SIMD instruction on ${SIMD_WIDTH} ${NUMBER_TYPE}s at a time")
add_executable(${TARGET_NAME} ./benchmark.cpp)
target_include_directories(${TARGET_NAME} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(${TARGET_NAME} mf-adr)
deal_ii_setup_target(${TARGET_NAME})
