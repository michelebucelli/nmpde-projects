#!/bin/bash

# PBS parameters for submission
# -j oe merges stdout and stderr logs
PBS_parameters="-q cpu -j oe"
number_repetition=10

# Handle plotting command directly
if [ "$1" == "--plot" ]; then
  apptainer exec "$(dirname "$0")/../dealii_paraview.sif" python3 "$(dirname "$0")/plot_scaling.py" "$2"
  exit $?
fi

# Check if we are running inside the compute node and parse parameters
is_node=false
first_job_id=""
desc=""
num_cpus=1
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
  elif [ "$1" == "--cpus" ]; then
    num_cpus="$2"
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
  OUTPUT_DIR="${first_job_id}_scalability"
  mkdir -p "${OUTPUT_DIR}"
  
  # Configure time format to output only elapsed seconds
  TIMEFORMAT="%R"
  
  # Execute the simulation inside the Apptainer container and measure execution time
  if [ "$num_cpus" -eq 1 ]; then
    # Run sequentially without mpirun
    elapsed=$( { time apptainer exec ../dealii_paraview.sif ./exercise-01 \
      $config_args \
      -o /dev/null \
      > "./${OUTPUT_DIR}/${PBS_JOBID}.out" \
      2> "./${OUTPUT_DIR}/${PBS_JOBID}.err"; } \
    2>&1 )
  else
    # Run with mpirun using the specified number of CPUs
    elapsed=$( { time mpirun -n $num_cpus apptainer exec ../dealii_paraview.sif ./exercise-01 \
      $config_args \
      -o /dev/null \
      > "./${OUTPUT_DIR}/${PBS_JOBID}.out" \
      2> "./${OUTPUT_DIR}/${PBS_JOBID}.err"; } \
    2>&1 )
  fi
  
  # File locking to prevent race conditions during concurrent writing
  results_file="${OUTPUT_DIR}/parallel_scalability.txt"
  lock_dir="${OUTPUT_DIR}/parallel_scalability.lock"
  
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
    if [ "$num_cpus" -eq 0 ]; then
      desc="0 CPUs (No MPI)"
    else
      desc="${num_cpus} CPUs"
    fi
  fi
  
  # Write the result row
  echo "$PBS_JOBID - $desc: ${elapsed} seconds" >> "$results_file"
  
else
  # ----------------------------------------------------
  # LOGIN NODE SUBMISSION
  # ----------------------------------------------------
  
  # Define CPU configurations to test (0 represents single process run without mpirun)
  cpu_configs=(1 2 4 8 16 32 56)
  procs_for_cpu=2

  SCRIPT=$(realpath $0)
  cd $(dirname $0)/..
  
  # Submit a separate PBS job for each CPU configuration
  first_job_id=""
  job_ids=""
  for num_cpus in "${cpu_configs[@]}"; do
    desc="${num_cpus}_MPI"
    for ((i=0; i<number_repetition; i++)); do
      qsub_cpus=$(( num_cpus / procs_for_cpu ))
      if [ "$num_cpus" -eq 1 ]; then
        qsub_cpus=1
        desc="No_MPI"
      fi
      
      if [ -z "$first_job_id" ]; then
        echo "Submitting first job with $num_cpus CPUs ($desc)"
        first_job_id=$(qsub -N "${desc}" $PBS_parameters -l ncpus=$qsub_cpus -- $SCRIPT --node-exec --cpus $num_cpus --desc "$desc")
        first_job_id=$(echo "$first_job_id" | tr -d '[:space:]')
        echo "First Job ID is: $first_job_id"
        job_ids="${first_job_id}"
      else
        echo "Submitting job with $num_cpus CPUs ($desc)"
        job_id=$(qsub -N "${desc}" $PBS_parameters -l ncpus=$qsub_cpus -- $SCRIPT --node-exec --first-job-id "$first_job_id" --cpus $num_cpus --desc "$desc")
        job_id=$(echo "$job_id" | tr -d '[:space:]')
        job_ids="${job_ids}:${job_id}"
      fi
    done
  done

  echo "Submitting dependent plotting job..."
  qsub -W depend=afterany:${job_ids} -N plot_cpu $PBS_parameters -l ncpus=1 -- $SCRIPT --plot "build/${first_job_id}_scalability"
fi
