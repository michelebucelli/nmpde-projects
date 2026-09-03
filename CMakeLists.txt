cmake_minimum_required(VERSION 3.13)
project(hybrid_matrix_free_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED)
find_package(deal.II REQUIRED
  HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR}
)
deal_ii_initialize_cached_variables()

include_directories(include)

add_executable(hybrid_matrix_free_solver src/main.cpp)
deal_ii_setup_target(hybrid_matrix_free_solver)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfloat-conversion -Wmissing-braces -Wnon-virtual-dtor")

enable_testing()

foreach(test test_problem_data test_operator test_mg_preconditioner test_parameters)
  add_executable(${test} tests/${test}.cpp)
  target_include_directories(${test} PRIVATE tests)
  deal_ii_setup_target(${test})
  add_test(NAME ${test} COMMAND ${test})
endforeach()
