#!/bin/bash

# PBS parameters for submission
# -l ncpus=28 requests 28 cores on the node
# -j oe merges stdout and stderr logs

PBS_parameters="-q cpu -l ncpus=28 -j oe"
number_repetition=10

# Handle plotting command directly
if [ "$1" == "--plot" ]; then
  apptainer exec "$(dirname "$0")/../dealii_paraview.sif" python3 "$(dirname "$0")/plot_performance.py" "$2"
  exit $?
fi

# Check if we are running inside the compute node and parse parameters
is_node=false
first_job_id=""
desc=""
config_args=""

# Process arguments
while [ $# -gt 0 ]; do
  if [ "$1" == "--node-exec" ]; then
    is_node=true
  elif [ "$1" == "--first-job-id" ]; then
    first_job_id="$2"
    shift
  elif [ "$1" == "--desc" ]; then
    desc="$2"
    shift
  else
    config_args="$config_args $1"
  fi
  shift
done

if [ "$is_node" = true ]; then
  # ----------------------------------------------------
  # COMPUTE NODE EXECUTION
  # ----------------------------------------------------
  
  # Move to the workspace directory where the job was submitted
  cd $PBS_O_WORKDIR/build
  
  # If no first_job_id was passed, default to our own PBS_JOBID
  if [ -z "$first_job_id" ]; then
    if [ -n "$PBS_JOBID" ]; then
      first_job_id="$PBS_JOBID"
    else
      first_job_id="local_$(date +%Y%m%d_%H%M%S)"
    fi
  fi
  
  # Create a dedicated directory for this job sweep
  OUTPUT_DIR="${first_job_id}_performance"
  mkdir -p "${OUTPUT_DIR}"
  
  # Configure time format to output only elapsed seconds
  TIMEFORMAT="%R"
  
  # Execute the simulation inside the Apptainer container and measure execution time
  elapsed=$( { time mpirun apptainer exec ../dealii_paraview.sif ./exercise-01 \
    $config_args \
    -o /dev/null \
    > "./${OUTPUT_DIR}/${PBS_JOBID}.out" \
    2> "./${OUTPUT_DIR}/${PBS_JOBID}.err"; } \
  2>&1 )
  
  # File locking to prevent race conditions during concurrent writing
  results_file="${OUTPUT_DIR}/solver_performance.txt"
  lock_dir="${OUTPUT_DIR}/solver_performance.lock"
  
  # Spin-lock using atomic directory creation
  while ! mkdir "$lock_dir" 2>/dev/null; do
    sleep 0.1
  done
  
  # Register trap to clean up the lock on exit (success or failure)
  trap 'rmdir "$lock_dir" 2>/dev/null' EXIT
  
  # Initialize the results file if it does not exist yet
  if [ ! -f "$results_file" ]; then
    echo "Configuration performance (JobID, Description, WallTime_Seconds):" > "$results_file"
    echo "------------------------------------------------------------------------" >> "$results_file"
  fi
  
  # If no desc was passed, fallback to the configuration arguments
  if [ -z "$desc" ]; then
    desc="$config_args"
  fi
  
  # Write the result row
  echo "$PBS_JOBID - $desc: ${elapsed} seconds" >> "$results_file"
  
else
  # ----------------------------------------------------
  # LOGIN NODE SUBMISSION
  # ----------------------------------------------------
  
  # Define configurations and their human-readable descriptions
  configurations=(\
    "-s direct" \
    "-s gmres -p ssor --ssor_omega 1.0" \
    "-s cg -p ssor --ssor_omega 1.0" \
    "-s cg -p ilu --ilu_fill 0" \
    "-s cg -p ilu --ilu_fill 2" \
    "-s cg -p ilu --ilu_fill 3" \
  )

  descriptions=(\
    "Direct Solver" \
    "GMRES + SSOR" \
    "CG + SSOR" \
    "CG + ILU(0)" \
    "CG + ILU(2)" \
    "CG + ILU(3)" \
  )

  SCRIPT=$(realpath $0)
  cd $(dirname $0)/..
  
  # Submit a separate PBS job for each configuration
  first_job_id=""
  job_ids=""
  for i in "${!configurations[@]}"; do
    config="${configurations[$i]}"
    desc="${descriptions[$i]}"
    for ((j=0; j<number_repetition; j++)); do
      
      if [ -z "$first_job_id" ]; then
        echo "Submitting first job for configuration: $config ($desc)"
        first_job_id=$(qsub -N $(echo $desc | tr ' ' '_' | tr -d '()') $PBS_parameters -- $SCRIPT $config --node-exec --desc "$desc")
        first_job_id=$(echo "$first_job_id" | tr -d '[:space:]')
        echo "First Job ID is: $first_job_id"
        job_ids="${first_job_id}"
      else
        echo "Submitting job for configuration: $config ($desc)"
        job_id=$(qsub -N $(echo $desc | tr ' ' '_' | tr -d '()' ) $PBS_parameters -- $SCRIPT $config --node-exec --first-job-id "$first_job_id" --desc "$desc")
        job_id=$(echo "$job_id" | tr -d '[:space:]')
        job_ids="${job_ids}:${job_id}"
      fi
    done
  done

  echo "Submitting dependent plotting job..."
  qsub -W depend=afterany:${job_ids} -N plot_solver $PBS_parameters -- $SCRIPT --plot "build/${first_job_id}_performance"
fi
