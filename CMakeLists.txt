cmake_minimum_required(VERSION 3.12.0)
project(project_title LANGUAGES CXX C)

include(cmake-common.cmake)
include(defaults.cmake)

add_custom_command(OUTPUT mesh
                  COMMAND mkdir mesh)

add_custom_command(OUTPUT models
                  COMMAND mkdir models)

add_custom_command(OUTPUT ${DEFAULT_BASE_MESH_FILE}.geo
                  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${DEFAULT_BASE_MESH_FILE}.geo ${DEFAULT_BASE_MESH_FILE}.geo
                  COMMAND touch ${DEFAULT_BASE_MESH_FILE}.geo
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${DEFAULT_BASE_MESH_FILE}.geo
                  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/mesh)

add_custom_command(OUTPUT ${DEFAULT_BASE_MESH_FILE}.msh
                  # COMMAND gmsh ${CMAKE_CURRENT_SOURCE_DIR}/mesh/mesh-square.geo -setnumber h 0.1 -save -o mesh/mesh-square.msh
                  COMMAND ${MESH_COMPILING_COMMAND} -setnumber h ${DEFAULT_MESH_SIZE} -o ${DEFAULT_BASE_MESH_FILE}.msh
                  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${DEFAULT_BASE_MESH_FILE}.geo)

add_custom_command(OUTPUT models/models.csv
                  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/models/models.csv models/models.csv 
                  COMMAND touch models/models.csv 
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/models/models.csv
                  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/models)

add_custom_target(mesh_file ALL
                  DEPENDS ${DEFAULT_BASE_MESH_FILE}.msh)


add_custom_target(models_file ALL
                  DEPENDS models/models.csv)

string(JOIN " " MESH_COMPILING_COMMAND_STR ${MESH_COMPILING_COMMAND} "")
add_compile_definitions(MESH_COMPILING_COMMAND="${MESH_COMPILING_COMMAND_STR}"s)
add_compile_definitions(DEFAULT_BASE_MESH_FILE="${DEFAULT_BASE_MESH_FILE}"s)
add_compile_definitions(DEFAULT_MESH_SIZE=${DEFAULT_MESH_SIZE})
add_compile_definitions(DEFAULT_MODEL="${DEFAULT_MODEL}"s)
add_compile_definitions(DEFAULT_DELTA_T=${DEFAULT_DELTA_T})
add_compile_definitions(DEFAULT_MAX_TIME=${DEFAULT_MAX_TIME})
add_compile_definitions(DEFAULT_THETA=${DEFAULT_THETA})
add_compile_definitions(DEFAULT_OUTPUT_FILE="${DEFAULT_OUTPUT_FILE}"s)
add_compile_definitions(DEFAULT_SOLVER="${DEFAULT_SOLVER}"s)
add_compile_definitions(DEFAULT_MAX_ITER=${DEFAULT_MAX_ITER})
add_compile_definitions(DEFAULT_TOLERANCE=${DEFAULT_TOLERANCE})
add_compile_definitions(DEFAULT_REDUCE_FACTOR=${DEFAULT_REDUCE_FACTOR})
add_compile_definitions(DEFAULT_PRECONDITIONER="${DEFAULT_PRECONDITIONER}"s)
add_compile_definitions(DEFAULT_ILU_FILL_IN=${DEFAULT_ILU_FILL_IN})
add_compile_definitions(DEFAULT_SSOR_OMEGA=${DEFAULT_SSOR_OMEGA})

add_executable(exercise-01 src/exercise-01.cpp src/args.cpp src/Current.cpp src/Solver.cpp src/Precondition.cpp)
deal_ii_setup_target(exercise-01)
