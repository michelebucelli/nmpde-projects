cmake_minimum_required(VERSION 3.12.0)
project(navier_stokes_time_dependent LANGUAGES CXX C)

include(common/cmake-common.cmake)

include_directories(include)

set(SOURCES
  src/NavierStokes.cpp
  src/main.cpp
)

add_executable(navier-stokes ${SOURCES})
deal_ii_setup_target(navier-stokes)
