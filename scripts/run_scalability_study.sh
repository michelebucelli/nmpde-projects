#!/usr/bin/env bash
set -euo pipefail

# Lightweight serial size study. It is not an MPI strong-scaling test.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
cd "${repo_root}"

output_root="${SCALABILITY_OUTPUT_DIR:-scalability-output}"
subdivisions="${SCALABILITY_NS:-12 24 48 96}"
final_time="${SCALABILITY_FINAL_TIME:-0.1}"
time_step="${SCALABILITY_TIME_STEP:-0.005}"

source "${script_dir}/common.sh"
load_course_modules

if [[ ! -x ./build/wave_fem ]]; then
  cmake -S . -B build
  cmake --build build -j"$(nproc)"
fi

mkdir -p "${output_root}"
printf 'N,dofs,cells,elapsed_seconds\n' > "${output_root}/scaling.csv"
for n in ${subdivisions}; do
  run_dir="${output_root}/N${n}"
  mkdir -p "${run_dir}"
  start=$(date +%s.%N)
  ./build/wave_fem "${n}" 1 "${final_time}" "${time_step}" 100 \
    "${run_dir}" free cg ssor > "${run_dir}/run.log"
  end=$(date +%s.%N)
  elapsed=$(awk -v start="${start}" -v end="${end}" 'BEGIN {print end-start}')
  read -r cells dofs < <(awk -F, 'NR==2 {print $4, $5}' "${run_dir}/parameters.csv")
  printf '%s,%s,%s,%s\n' "${n}" "${dofs}" "${cells}" "${elapsed}" >> "${output_root}/scaling.csv"
done

echo "Wrote ${output_root}/scaling.csv"
echo "This is a serial problem-size study, not MPI strong or weak scaling."
