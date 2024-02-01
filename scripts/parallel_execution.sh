#!/bin/bash


PRECONDITIONERS=("identity" "jacobi" "ssor" "amg")
MESH="mesh-cube-40.msh"


for n in {1..4}
do
  for PRECONDITIONER in "${PRECONDITIONERS[@]}"
  do
    for ((c=1; c<=1; c++))
    do
      mpiexec -n "$n" ../build/Project_parallel -m "$MESH" -P "$PRECONDITIONER" -p 5
    done
  done
done


# python3 cond_number_vs_p_value.py
# python3 iterations_vs_p_value.py