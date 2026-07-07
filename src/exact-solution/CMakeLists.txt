cmake_minimum_required(VERSION 3.12.0)
project(space-time-adaptivity LANGUAGES CXX C)

include(common/cmake-common.cmake)

# Space-Time Adaptivity
add_executable(space-time-adaptivity 
    src/space-time-adaptivity/space-time-adaptivity.cpp 
    src/space-time-adaptivity/Heat.cpp
)
deal_ii_setup_target(space-time-adaptivity)

# Non-Adaptive
add_executable(non-adaptive 
    src/non-adaptive/non-adaptive.cpp 
    src/non-adaptive/Heat.cpp
)
deal_ii_setup_target(non-adaptive)

# Exact-Test
add_executable(exact-test 
    src/exact-test/exact-test.cpp 
    src/exact-test/Heat.cpp
)
deal_ii_setup_target(exact-test)

enable_testing()

add_test(NAME space-time-adaptivity
         COMMAND space-time-adaptivity)

add_test(NAME non-adaptive
         COMMAND non-adaptive)

add_test(NAME exact-test
         COMMAND exact-test)

add_test(NAME model_validation
         COMMAND space-time-adaptivity)