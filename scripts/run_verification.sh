#!/usr/bin/env bash
set -euo pipefail

# Simple verification runner: configure, build, run ctest
BUILD_DIR=${1:-build}

echo "Configuring (build dir: $BUILD_DIR)"
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DDEAL_II_ALLOW_TESTS=ON

echo "Building"
cmake --build "$BUILD_DIR" --parallel

echo "Running ctest"
cd "$BUILD_DIR"
ctest --output-on-failure --parallel $(( $(nproc) - 1 ))

echo "Verification complete"
