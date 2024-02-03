#!/bin/bash


PRECONDITIONERS=("identity" "ssor" ) # add to the vector new preconditioners

MESH="mesh-cube-40.msh"

for PRECONDITIONER in "${PRECONDITIONERS[@]}"
do
  for P_VALUE in {0..8}
  do
    for ((c=1; c<=1; c++))
    do
      mpiexec -n 4 ../build/Project_parallel -m "$MESH"  -P "$PRECONDITIONER" -p "$P_VALUE"
    done
  done
done

# python3 cond_number_vs_p_value.py
# python3 iterations_vs_p_value.py