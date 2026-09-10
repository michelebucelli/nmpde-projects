cmake_minimum_required(VERSION 3.12.0)
project(lab05 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

# Create a static library from ADR.cpp
add_library(adr_lib STATIC src/TimeADR.cpp)
deal_ii_setup_target(adr_lib)

add_executable(main src/main.cpp)
target_link_libraries(main adr_lib)
deal_ii_setup_target(main)

add_executable(convergence src/convergence.cpp)
target_link_libraries(convergence adr_lib)
deal_ii_setup_target(convergence)