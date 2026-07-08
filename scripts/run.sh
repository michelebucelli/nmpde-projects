#!/bin/bash

# PBS parameters for interactive job
PBS_parameters="-q cpu -l ncpus=28 -j oe"

# Check arguments
is_node=false
is_init=false
for arg in "$@"; do
  if [ "$arg" == "--node-exec" ]; then
    is_node=true
  elif [ "$arg" == "--init" ]; then
    is_init=true
  fi
done

if [ "$is_node" = true ]; then
  # ----------------------------------------------------
  # COMPUTE NODE EXECUTION
  # ----------------------------------------------------
  
  # Get project root (one level out from script directory)
  SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
  PROJECT_ROOT=$(dirname "$SCRIPT_DIR")
  cd "$PROJECT_ROOT"
  
  if [ ! -f "dealii_paraview.sif" ]; then
    echo "=== Starting Apptainer Image Build ==="
    # Build the Apptainer image
    apptainer build --fakeroot dealii_paraview.sif dealii_paraview.def
    if [ $? -ne 0 ]; then
      echo "Error: Apptainer build failed!"
      exit 1
    fi
    echo "=== Apptainer Image Build Completed ==="
  else
    echo "=== Apptainer Image dealii_paraview.sif already exists. Skipping build. ==="
  fi

  echo "=== Starting Project Compilation ==="
  mkdir -p build
  cd build
  
  # Run cmake and make inside the new apptainer image
  apptainer exec ../dealii_paraview.sif bash -c "cmake .. && make -j\$(nproc)"
  if [ $? -ne 0 ]; then
    echo "Error: Compilation failed!"
    exit 1
  fi
  echo "=== Project Compilation Completed ==="

else
  # ----------------------------------------------------
  # LOGIN NODE SUBMISSION
  # ----------------------------------------------------
  SCRIPT=$(realpath "$0")
  PROJECT_ROOT=$(dirname "$(dirname "$SCRIPT")")
  cd "$PROJECT_ROOT"
  
  echo "Requesting interactive PBS job on compute node..."
  # Submit an interactive job (-I) passing the script command after --
  qsub -I $PBS_parameters -- "$SCRIPT" --node-exec
  
  # If --init flag was specified, do not execute the remaining sweep scripts
  if [ "$is_init" = true ]; then
    echo "=== Initialization finished. Skipping sweep scripts execution. ==="
    exit 0
  fi

  echo "=== Triggering all other sweep scripts ==="
  for script in scripts/*.sh; do
    # Check if the file exists and is not this script
    if [ -f "$script" ] && [ "$(basename "$script")" != "$(basename "$0")" ]; then
      echo "Running script: $script"
      chmod +x "$script" 2>/dev/null
      ./"$script"
    fi
  done

  echo "=== Build and Sweep Initialization Finished ==="
fi
