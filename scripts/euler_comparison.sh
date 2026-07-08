#!/bin/bash

# PBS parameters for submission
# -l ncpus=28 requests 28 cores on the node
# -j oe merges stdout and stderr logs
PBS_parameters="-q cpu -l ncpus=28 -j oe"

# Handle plotting command directly
if [ "$1" == "--plot" ]; then
  pvpython "$(dirname "$0")/plot_euler_comparison.py" "${@:2}"
  exit $?
fi

# Check if we are running inside the compute node and parse parameters
is_node=false
method_flag=""
method_name=""

# Process arguments
while [ $# -gt 0 ]; do
  if [ "$1" == "--node-exec" ]; then
    is_node=true
  elif [ "$1" == "--explicit_euler" ] || [ "$1" == "--crank_nicolson" ] || [ "$1" == "--implicit_euler" ]; then
    method_flag="$1"
    method_name="${1#--}" # Strip leading --
  fi
  shift
done

if [ "$is_node" = true ]; then
  # ----------------------------------------------------
  # COMPUTE NODE EXECUTION
  # ----------------------------------------------------
  
  if [ -z "$method_flag" ]; then
    echo "Error: No time integration method specified."
    exit 1
  fi

  # Move to the workspace directory where the job was submitted
  cd $PBS_O_WORKDIR/build
  
  # Determine job ID
  job_id="$PBS_JOBID"
  if [ -z "$job_id" ]; then
    job_id="local_$(date +%Y%m%d_%H%M%S)"
  fi
  
  # Define and create folder-specific directory inside build/
  OUTPUT_DIR="${method_name}_${job_id}"
  mkdir -p "$OUTPUT_DIR"
  
  echo "Running simulation for integration method: $method_name"
  
  # Run the simulation inside the Apptainer container
  # Save the activation time output inside the method_jobid directory
  mpirun apptainer exec ../dealii_paraview.sif ./exercise-01 \
    "$method_flag" \
    -o "${OUTPUT_DIR}/activation_time"
    
  echo "Finished simulation for: $method_name"
  
else
  # ----------------------------------------------------
  # LOGIN NODE SUBMISSION
  # ----------------------------------------------------
  
  methods=(\
    "explicit_euler" \
    "crank_nicolson" \
    "implicit_euler" \
  )
  
  SCRIPT=$(realpath $0)
  cd $(dirname $0)/..
  pwd
  echo $SCRIPT
  
  # Submit a separate PBS job for each time integration method
  job_ids=""
  for method in "${methods[@]}"; do
    echo "Submitting job for method: $method"
    job_id=$(qsub -N $method $PBS_parameters -- $SCRIPT "--${method}" --node-exec)
    job_id=$(echo "$job_id" | tr -d '[:space:]')
    if [ -z "$job_ids" ]; then
      job_ids="${job_id}"
    else
      job_ids="${job_ids}:${job_id}"
    fi
  done

  echo "Submitting dependent plotting job..."
  qsub -W depend=afterany:${job_ids} -N plot_euler $PBS_parameters -- $SCRIPT --plot
fi
