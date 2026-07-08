#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   TEST_EXEC=build/space-time-adaptivity ./scripts/run_scalability.sh build
# space-time-adaptivity is MPI-parallel (not thread-based), so scaling is
# done via `mpirun -np N`, not a --threads flag.
# If TEST_EXEC is not set, the script will try to run ctest -R scalability

BUILD_DIR=${1:-build}
TEST_EXEC=${TEST_EXEC:-}
OUT_FILE=${OUT_FILE:-scalability_results.csv}
PROCESS_COUNTS=${PROCESS_COUNTS:-"1 2 4 8 16"}

if [ -n "$TEST_EXEC" ]; then
  echo "processes,seconds,cmd" > "$OUT_FILE"
  for np in $PROCESS_COUNTS; do
    echo "Running $TEST_EXEC with $np MPI process(es)"
    START=$(date +%s.%N)
    mpirun -np "$np" "$TEST_EXEC"
    END=$(date +%s.%N)
    DIFF=$(echo "$END - $START" | bc)
    echo "$np,$DIFF,mpirun -np $np $TEST_EXEC" >> "$OUT_FILE"
  done
  echo "Results written to $OUT_FILE"
  cat "$OUT_FILE"
else
  echo "No TEST_EXEC provided, falling back to running ctest -R scalability (in $BUILD_DIR)"
  cd "$BUILD_DIR"
  ctest -R scalability --output-on-failure
fi
