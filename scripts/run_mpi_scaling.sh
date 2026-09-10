#!/usr/bin/env bash
set -euo pipefail

# Measure the same MPI problem with several process counts.
source "$(dirname "$0")/common.sh"
load_course_modules

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
output_root="${MPI_SCALING_OUTPUT_DIR:-${repo_root}/output/mpi-scaling}"
mpi_procs="${MPI_PROCS:-1 2 4}"
n="${MPI_N:-64}"
final_time="${MPI_FINAL_TIME:-0.05}"
dt="${MPI_TIME_STEP:-0.005}"
case_name="${MPI_CASE:-free}"
steps="$(awk -v T="${final_time}" -v dt="${dt}" 'BEGIN { printf "%d", T/dt }')"

if [[ ! -x "${repo_root}/build/wave_fem" || ! -x "${repo_root}/build/wave_fem_mpi" ]]; then
  echo "Missing build/wave_fem or build/wave_fem_mpi; build the project first." >&2
  exit 1
fi

mkdir -p "${output_root}"
summary="${output_root}/scaling.csv"
printf 'implementation,n_ranks,n_subdivisions,dofs,cells,setup_seconds,time_loop_seconds,postprocess_seconds,elapsed_seconds\n' > "${summary}"

serial_dir="${output_root}/serial"
rm -rf "${serial_dir}"
echo "Running serial baseline..."
mkdir -p "${serial_dir}"
serial_start_ns="$(date +%s%N)"
if "${repo_root}/build/wave_fem" "${n}" 1 "${final_time}" "${dt}" "$((steps + 1))" \
    "${serial_dir}" "${case_name}" >/dev/null; then
  serial_end_ns="$(date +%s%N)"
  serial_elapsed="$(awk -v start="${serial_start_ns}" -v end="${serial_end_ns}" \
    'BEGIN { printf "%.9f", (end-start)/1000000000 }')"
else
  echo "Serial baseline failed." >&2
  exit 1
fi
printf 'serial,1,%s,%s,%s,0,0,0,%s\n' "${n}" "$(( (n + 1) * (n + 1) ))" \
  "$(( n * n ))" "${serial_elapsed}" >> "${summary}"

for p in ${mpi_procs}; do
  run_dir="${output_root}/np${p}"
  rm -rf "${run_dir}"
  echo "Running ${p} MPI process(es)..."
  if [[ "$(id -u)" -eq 0 ]]; then
    mpirun --allow-run-as-root -np "${p}" "${repo_root}/build/wave_fem_mpi" \
      "${n}" "${final_time}" "${dt}" "${run_dir}" "${case_name}"
  else
    mpirun -np "${p}" "${repo_root}/build/wave_fem_mpi" \
      "${n}" "${final_time}" "${dt}" "${run_dir}" "${case_name}"
  fi
  printf 'mpi,' >> "${summary}"
  tail -n 1 "${run_dir}/timing.csv" >> "${summary}"
done

echo "MPI scaling results: ${summary}"
