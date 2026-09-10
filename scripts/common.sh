#!/usr/bin/env bash

# Shared helpers for the course Docker/HPC environment.
# This file is meant to be sourced by the other shell scripts.

load_course_modules() {
  if [[ -f /u/sw/etc/profile ]]; then
    # The Docker image used for the course does not expose `module` until this
    # profile is sourced.
    export LOGNAME="${LOGNAME:-$(whoami)}"
    source /u/sw/etc/profile
  fi

  if command -v module >/dev/null 2>&1; then
    module load gcc-glibc/11.2.0 dealii/9.5.1
  fi
}

plot_diagnostics_if_possible() {
  local csv_file="$1"

  if ! command -v python >/dev/null 2>&1; then
    echo "python not available; skipping PNG plots."
    return
  fi

  if python -c "import matplotlib" >/dev/null 2>&1; then
    python scripts/plot_diagnostics.py "${csv_file}"
  else
    echo "matplotlib not available; skipping PNG plots."
    echo "Install if needed with: python -m pip install --user matplotlib"
  fi
}
