cmake_minimum_required(VERSION 3.12.0)
project(project_title LANGUAGES CXX C)

include(common/cmake-common.cmake)

add_executable(H_Heat src/homogeneous/main.cpp src/homogeneous/H_Heat.cpp src/timing/Timing.cpp)
deal_ii_setup_target(H_Heat)

add_executable(STA_Heat src/adaptive/main.cpp src/adaptive/STA_Heat.cpp src/homogeneous/H_Heat.cpp src/timing/Timing.cpp)
deal_ii_setup_target(STA_Heat)
