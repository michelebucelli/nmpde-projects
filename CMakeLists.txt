cmake_minimum_required(VERSION 3.12.0)
project(project_title LANGUAGES CXX C)

include(./common/cmake-common.cmake)

set(CMAKE_BUILD_TYPE Release)  # Ensure Release mode
set(CMAKE_CXX_FLAGS_RELEASE "-O3")  # Explicitly set -O3 for Release mode
set(CMAKE_C_FLAGS_RELEASE "-O3")

add_executable(main src/main.cpp src/UncoupledNavierStokes.cpp src/MonolithicNavierStokes.cpp src/SteadyNavierStokes.cpp src/ConfigReader.cpp)
deal_ii_setup_target(main)
