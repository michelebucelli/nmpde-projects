#!/bin/bash

# PBS parameters for submission
# -l ncpus=28 requests 28 cores on the node
# -j oe merges stdout and stderr logs
PBS_parameters="-q cpu -l ncpus=28 -j oe"

# Check if we are running inside the compute node and parse parameters
is_node=false
first_job_id=""
model=""
run_plot=false
plot_dir=""

# Process arguments
while [ $# -gt 0 ]; do
  if [ "$1" == "--node-exec" ]; then
    is_node=true
  elif [ "$1" == "--first-job-id" ]; then
    first_job_id="$2"
    shift
  elif [ "$1" == "-m" ]; then
    model="$2"
    shift
  elif [ "$1" == "--plot" ]; then
    run_plot=true
    if [ $# -gt 1 ] && [[ "$2" != -* ]]; then
      plot_dir="$2"
      shift
    fi
  fi
  shift
done

if [ "$run_plot" = true ]; then
  SCRIPT_DIR=$(dirname $(realpath $0))
  echo "Running model clip rendering inside Apptainer..."
  apptainer exec "$SCRIPT_DIR/../dealii_paraview.sif" pvpython "$SCRIPT_DIR/plot_model_activation.py" "$plot_dir"
  exit $?
fi

if [ "$is_node" = true ]; then
  # ----------------------------------------------------
  # COMPUTE NODE EXECUTION
  # ----------------------------------------------------
  
  if [ -z "$model" ]; then
    echo "Error: No model specified."
    exit 1
  fi

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
  
  # Define and create model-specific directory nested inside common parent inside build/
  OUTPUT_DIR="${first_job_id}_models/${model}_${job_id}"
  mkdir -p "$OUTPUT_DIR"
  
  echo "Running simulation for model: $model"
  
  # Run the simulation inside the Apptainer container
  # Save the activation time output and redirect stdout/stderr inside the directory
  mpirun apptainer exec ../dealii_paraview.sif ./exercise-01 \
    -m "$model" \
    -o "${OUTPUT_DIR}/${model}_activation_time" \
    > "${OUTPUT_DIR}/${job_id}.out" \
    2> "${OUTPUT_DIR}/${job_id}.err"
    
  echo "Finished model: $model"
  
else
  # ----------------------------------------------------
  # LOGIN NODE SUBMISSION
  # ----------------------------------------------------
  
  MODELS_FILE="models/models.csv"
  if [ ! -f "$MODELS_FILE" ]; then
    echo "Error: $MODELS_FILE not found."
    exit 1
  fi
  
  # Read model names from the CSV file (first column, skipping the header)
  models=()
  first_line=true
  while IFS=, read -r col1 remainder || [ -n "$col1" ]; do
    if [ "$first_line" = true ]; then
      first_line=false
      continue
    fi
    # Trim whitespace and quotes
    model_name=$(echo "$col1" | tr -d '[:space:]"' )
    if [ -n "$model_name" ]; then
      models+=("$model_name")
    fi
  done < "$MODELS_FILE"
  
  SCRIPT=$(realpath $0)
  cd $(dirname $0)/..
  
  # Submit a separate PBS job for each model
  first_job_id=""
  job_ids=""
  for m in "${models[@]}"; do
    if [ -z "$first_job_id" ]; then
      echo "Submitting first job for model: $m"
      first_job_id=$(qsub -N ${m} $PBS_parameters -- $SCRIPT -m "$m" --node-exec)
      first_job_id=$(echo "$first_job_id" | tr -d '[:space:]')
      echo "First Job ID is: $first_job_id"
      job_ids="${first_job_id}"
    else
      echo "Submitting job for model: $m"
      job_id=$(qsub -N ${m} $PBS_parameters -- $SCRIPT -m "$m" --node-exec --first-job-id "$first_job_id")
      job_id=$(echo "$job_id" | tr -d '[:space:]')
      job_ids="${job_ids}:${job_id}"
    fi
  done

  echo "Submitting dependent plotting job..."
  qsub -W depend=afterany:${job_ids} -N plot_model $PBS_parameters -- $SCRIPT --plot "build/${first_job_id}_models"
fi
