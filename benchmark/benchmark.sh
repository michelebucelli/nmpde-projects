DIMS=(2)
DEGREES=(1 2 3)
MPI_PROCESSES=(2 4 8)
REAL_NUMBERS=("double" "float")
SIMD_WIDTHS=(2 4 8)

for DIM in "${DIMS[@]}"; do
    for DEG in "${DEGREES[@]}"; do
        for RN in "${REAL_NUMBERS[@]}"; do
            for SW in "${SIMD_WIDTHS[@]}"; do
                if [[ "$RN" == "double" && "$SW" == 8 ]] || [[ "$RN" == "float" && "$SW" == 2 ]]; then
                    continue 
                fi
                CURRENT_BUILD_DIR="./builds/build_d${DIM}_p${DEG}_w${SW}_${RN}"

                cmake -B ${CURRENT_BUILD_DIR} -DFE_DEGREE=${DEG} -DPROBLEM_DIMENSIONS=${DIM} -DNUMBER_TYPE=${RN} -DSIMD_WIDTH=${SW} . \
                    || { echo "--> CMake configuration failed. Skipping."; continue; }
                cmake --build ${CURRENT_BUILD_DIR} -j$(nproc) \
                    || { echo "--> Build failed. Skipping."; continue; }

                for p in "${MPI_PROCESSES[@]}"; do                
                    echo "Execution (${p} processes) Matrix Free ADR with ${DIM}D mesh and FE of degree ${DEG}, using ${SW} ${RN} SIMD instructions"
                    mpirun -n ${p} ${CURRENT_BUILD_DIR}/bex_d${DIM}_p${DEG}_w${SW}_${RN}
                done
            done
        done
    done
done
