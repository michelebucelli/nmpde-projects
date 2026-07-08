cmake_minimum_required(VERSION 3.12.0)
project(NavierStokes LANGUAGES CXX C)

include(common/cmake-common.cmake)

add_executable(main2D src/main2D.cpp src/NavierStokes2D.cpp)
deal_ii_setup_target(main2D)

add_executable(main3D src/main3D.cpp src/NavierStokes3D.cpp)
deal_ii_setup_target(main3D)