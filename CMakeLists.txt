cmake_minimum_required(VERSION 3.12.0)
project(space-time-adaptivity LANGUAGES CXX C)

include(common/cmake-common.cmake)

add_executable(space-time-adaptivity src/space-time-adaptivity.cpp src/Heat.cpp)
deal_ii_setup_target(space-time-adaptivity)

add_executable(manufactured-solution-test src/manufactured-solution-test.cpp src/Heat.cpp)
deal_ii_setup_target(manufactured-solution-test)

enable_testing()

add_test(NAME space-time-adaptivity
         COMMAND space-time-adaptivity)
set_tests_properties(space-time-adaptivity PROPERTIES FIXTURES_SETUP space_time_output)

add_test(NAME manufactured-solution-test
         COMMAND manufactured-solution-test)

add_test(NAME model_validation
         COMMAND python3 ${CMAKE_SOURCE_DIR}/scripts/validate_model.py
                 --results ${CMAKE_CURRENT_BINARY_DIR}/mass_balance-n12-adaptive.csv
                 --tolerance 1e-2)
set_tests_properties(model_validation PROPERTIES FIXTURES_REQUIRED space_time_output)

set_tests_properties(manufactured-solution-test PROPERTIES DEPENDS model_validation)
