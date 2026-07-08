#!/bin/bash

# PBS parameters for submission
# -l ncpus=28 requests 28 cores on the node
# -j oe merges stdout and stderr logs
PBS_parameters="-q cpu -l ncpus=28 -j oe"

# Check if we are running inside the compute node and parse parameters
is_node=false
first_job_id=""
mesh_size="0.1"
delta_t="0.005"
desc=""
run_plot=false
plot_dir=""

# Process arguments
while [ $# -gt 0 ]; do
  if [ "$1" == "--node-exec" ]; then
    is_node=true
  elif [ "$1" == "--first-job-id" ]; then
    first_job_id="$2"
    shift
  elif [ "$1" == "--mesh_size" ]; then
    mesh_size="$2"
    shift
  elif [ "$1" == "--delta_t" ]; then
    delta_t="$2"
    shift
  elif [ "$1" == "--desc" ]; then
    desc="$2"
    shift
  elif [ "$1" == "--plot" ]; then
    run_plot=true
    if [ $# -gt 1 ] && [[ "$2" != --* ]]; then
      plot_dir="$2"
      shift
    fi
  fi
  shift
done

if [ "$run_plot" = true ]; then
  # Determine script directory
  SCRIPT_DIR=$(dirname $(realpath $0))
  echo "Running 3D plot generation inside Apptainer..."
  apptainer exec "$SCRIPT_DIR/../dealii_paraview.sif" pvpython "$SCRIPT_DIR/plot_convergence_3d.py" "$plot_dir"
  exit $?
fi

if [ "$is_node" = true ]; then
  # ----------------------------------------------------
  # COMPUTE NODE EXECUTION
  # ----------------------------------------------------
  
  # Move to the workspace directory where the job was submitted
  cd $PBS_O_WORKDIR/build
  
  # Determine job ID
  job_id="$PBS_JOBID"
  if [ -z "$job_id" ]; then
    job_id="local_$(date +%Y%m%d_%H%M%S)"
  fi
  
  # Default first_job_id to our own if not passed (e.g. first job run)
  if [ -z "$first_job_id" ]; then
    first_job_id="$job_id"
  fi
  
  # Generate description if not passed
  if [ -z "$desc" ]; then
    h_str=$(echo "$mesh_size" | tr '.' '_')
    dt_str=$(echo "$delta_t" | tr '.' '_')
    desc="h${h_str}_dt${dt_str}"
  fi
  
  # Create a subdirectory for this run nested under the common parent sweep directory
  OUTPUT_DIR="${first_job_id}_convergence/${desc}_${job_id}"
  mkdir -p "$OUTPUT_DIR"
  
  echo "Running simulation for: $desc (h=${mesh_size}mm, dt=${delta_t}ms)"
  
  # Run the simulation inside the Apptainer container
  # Redirect stdout and stderr to logs inside the subdirectory
  mpirun apptainer exec ../dealii_paraview.sif ./exercise-01 \
    --mesh_size "$mesh_size" \
    --delta_t "$delta_t" \
    -o "${OUTPUT_DIR}/activation_time" \
    > "${OUTPUT_DIR}/${job_id}.out" \
    2> "${OUTPUT_DIR}/${job_id}.err"
    
  echo "Finished simulation for: $desc"
  
else
  # ----------------------------------------------------
  # LOGIN NODE SUBMISSION
  # ----------------------------------------------------
  
  mesh_sizes=(0.1 0.2 0.5)
  delta_ts=(0.005 0.01 0.05)
  
  SCRIPT=$(realpath $0)
  cd $(dirname $0)/..
  pwd
  echo $SCRIPT
  
  # Submit a separate PBS job for each combination of mesh_size and delta_t
  first_job_id=""
  job_ids=""
  for h in "${mesh_sizes[@]}"; do
    for dt in "${delta_ts[@]}"; do
      # Substitute dots with underscores for job name safety
      h_str=$(echo "$h" | tr '.' '_')
      dt_str=$(echo "$dt" | tr '.' '_')
      desc="h${h_str}__dt${dt_str}"
      
      if [ -z "$first_job_id" ]; then
        echo "Submitting first job for combination: h: ${h}, dt: ${dt}"
        first_job_id=$(qsub -N ${desc} $PBS_parameters -- $SCRIPT --mesh_size "$h" --delta_t "$dt" --node-exec --desc "$desc")
        first_job_id=$(echo "$first_job_id" | tr -d '[:space:]')
        echo "First Job ID is: $first_job_id"
        job_ids="${first_job_id}"
      else
        echo "Submitting job for combination: h: ${h}, dt: ${dt}"
        job_id=$(qsub -N ${desc} $PBS_parameters -- $SCRIPT --mesh_size "$h" --delta_t "$dt" --node-exec --first-job-id "$first_job_id" --desc "$desc")
        job_id=$(echo "$job_id" | tr -d '[:space:]')
        job_ids="${job_ids}:${job_id}"
      fi
    done
  done

  echo "Submitting dependent plotting job..."
  qsub -W depend=afterany:${job_ids} -N plot_conv $PBS_parameters -- $SCRIPT --plot "build/${first_job_id}_convergence"
fi
