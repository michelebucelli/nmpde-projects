#!/usr/bin/env bash
set -euo pipefail

# run_model_validation.sh [build_dir]
VENV_DIR=${VENV_DIR:-venv}
BUILD_DIR=${1:-build}
RESULTS_FILE=${RESULTS_FILE:-model_results.csv}
TEST_CMD=${TEST_CMD:-}
PY_VALIDATOR=${PY_VALIDATOR:-scripts/validate_model.py}
TOL=${TOL:-1e-2}

echo "Build dir: $BUILD_DIR"

# Set up Python virtual environment and install dependencies
if [ ! -d "$VENV_DIR" ]; then
  echo "Creating Python virtual environment in $VENV_DIR"
  python3 -m venv "$VENV_DIR"
fi

echo "Installing/updating validation dependencies..."
# Run pip in an isolated subshell to ignore external PYTHONPATH
(
  unset PYTHONPATH
  . "$VENV_DIR/bin/activate"
  pip install --upgrade pip
  pip install -r scripts/requirements.txt
)

if [ -n "$TEST_CMD" ]; then
  echo "Running provided test command: $TEST_CMD"
  eval "$TEST_CMD" > "$RESULTS_FILE"
  echo "Wrote results to $RESULTS_FILE"
else
  echo "No TEST_CMD provided, attempting to run ctest -R space-time-adaptivity"
  if ! (cd "$BUILD_DIR" && ctest -R '^space-time-adaptivity$' --output-on-failure); then
    echo "ctest command failed. Aborting validation." >&2
    exit 1
  fi

  RESULTS_FILE="$BUILD_DIR/mass_balance-n12-adaptive.csv"
  if [ ! -e "$RESULTS_FILE" ]; then
    echo "No mass_balance-n12-adaptive.csv found at $RESULTS_FILE; exiting with error"
    exit 2
  fi
  echo "Found mass-balance results: $RESULTS_FILE"
fi
echo "Running Python validator ($PY_VALIDATOR) with tolerance $TOL"

# Run the validator in an isolated subshell using the virtual environment
(
  unset PYTHONPATH
  . "$VENV_DIR/bin/activate"
  python3 "$PY_VALIDATOR" --results "$RESULTS_FILE" --tolerance "$TOL"
)

echo "Model validation completed"
