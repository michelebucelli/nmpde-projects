cmake_minimum_required(VERSION 3.12.0)
project(lab08 LANGUAGES CXX C)

include(../common/cmake-common.cmake)



add_executable(FK_solver src/main.cpp src/FK_solver.cpp)
deal_ii_setup_target(FK_solver)
