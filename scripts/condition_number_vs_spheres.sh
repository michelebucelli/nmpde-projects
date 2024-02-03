#!/bin/bash


PRECONDITIONERS=("identity" "jacobi" "ssor" "amg")
MESH="mesh-cube-40.msh"


for n in {1..7}
do
  for PRECONDITIONER in "${PRECONDITIONERS[@]}"
  do
    for ((c=1; c<=1; c++))
    do
      mpiexec -n 4 ../build/Project_parallel -m "$MESH" -P "$PRECONDITIONER" -p 2 -n "$n"
    done
  done
done


# python3 cond_number_vs_p_value.py
# python3 iterations_vs_p_value.py