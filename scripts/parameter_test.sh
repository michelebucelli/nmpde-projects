#!/bin/bash


SWEEPS=("1" "5" "10" "15" "20")
PRECONDITIONER="ssor" # add to the vector new preconditioners
MESH="mesh-cube-20.msh"


for SWEEP in "${SWEEPS[@]}"
do
  for P_VALUE in {0..8}
  do
    for ((c=1; c<=1; c++))
    do
      mpiexec -n 4 ../build/Project_parallel -m "$MESH" -P "$PRECONDITIONER" -p "$P_VALUE" -o "$SWEEP"
    done
  done
done


# python3 cond_number_vs_p_value.py
# python3 iterations_vs_p_value.py