#!/usr/bin/env bash
set -euo pipefail

# Small reproducibility helper for the report.
#
# It assumes that the project has already been built and that ./build/wave_fem
# exists. Each run writes to a different output folder so results are not
# overwritten.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

source "${script_dir}/common.sh"

usage() {
  cat <<EOF
Usage:
  bash scripts/run_validation_suite.sh

Optional environment variables:
  VALIDATION_OUTPUT_DIR   Output root folder. Default: output/validation/suite
  SPACE_NS                Space-separated mesh subdivisions. Default: "12 24 48"
  SPACE_DEGREE            FE degree for space checks. Default: 1
  SPACE_FINAL_TIME        Final time for space checks. Default: 1.0
  SPACE_TIME_STEP         Time step for space checks. Default: 0.0025
  SPACE_OUTPUT_FREQUENCY  Output frequency for space checks. Default: 100
  TIME_N_SUBDIVISIONS     Mesh subdivisions for time checks. Default: 32
  TIME_DEGREE             FE degree for time checks. Default: 1
  TIME_FINAL_TIME         Final time for time checks. Default: 1.0
  TIME_DTS                Space-separated time steps. Default: "0.02 0.01 0.005"
  TIME_OUTPUT_FREQUENCY   Output frequency for time checks. Default: 100
  VALIDATION_CASES        Space-separated cases. Default: "free forced"
  MAKE_PLOTS              yes/no. Default: yes

Example:
  SPACE_NS="16 32" TIME_DTS="0.01 0.005" VALIDATION_CASES="free forced" bash scripts/run_validation_suite.sh
EOF
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

cd "${repo_root}"

executable="./build/wave_fem"
validation_output_dir="${VALIDATION_OUTPUT_DIR:-output/validation/suite}"
space_ns="${SPACE_NS:-12 24 48}"
space_degree="${SPACE_DEGREE:-1}"
space_final_time="${SPACE_FINAL_TIME:-1.0}"
space_time_step="${SPACE_TIME_STEP:-0.0025}"
space_output_frequency="${SPACE_OUTPUT_FREQUENCY:-100}"
time_n_subdivisions="${TIME_N_SUBDIVISIONS:-32}"
time_degree="${TIME_DEGREE:-1}"
time_final_time="${TIME_FINAL_TIME:-1.0}"
time_dts="${TIME_DTS:-0.02 0.01 0.005}"
time_output_frequency="${TIME_OUTPUT_FREQUENCY:-100}"
validation_cases="${VALIDATION_CASES:-free forced}"
make_plots="${MAKE_PLOTS:-yes}"

load_course_modules

if [[ ! -x "${executable}" ]]; then
  echo "Missing ${executable}. Build first with:"
  echo "  cmake -S . -B build"
  echo "  cmake --build build -j\$(nproc)"
  exit 1
fi

mkdir -p "${validation_output_dir}"

for validation_case in ${validation_cases}; do
  echo "Running validation case: ${validation_case}"

  echo "Running space refinement checks"
  for n in ${space_ns}; do
    output_dir="${validation_output_dir}/${validation_case}/space_N${n}"
    "${executable}" \
      "${n}" \
      "${space_degree}" \
      "${space_final_time}" \
      "${space_time_step}" \
      "${space_output_frequency}" \
      "${output_dir}" \
      "${validation_case}"

    if [[ "${make_plots}" == "yes" ]]; then
      plot_diagnostics_if_possible "${output_dir}/diagnostics.csv"
    fi
  done

  echo "Running time-step refinement checks"
  for dt in ${time_dts}; do
    safe_dt="${dt/./p}"
    output_dir="${validation_output_dir}/${validation_case}/time_dt${safe_dt}"
    "${executable}" \
      "${time_n_subdivisions}" \
      "${time_degree}" \
      "${time_final_time}" \
      "${dt}" \
      "${time_output_frequency}" \
      "${output_dir}" \
      "${validation_case}"

    if [[ "${make_plots}" == "yes" ]]; then
      plot_diagnostics_if_possible "${output_dir}/diagnostics.csv"
    fi
  done
done

echo "Validation suite completed. Inspect:"
find "${validation_output_dir}" -name diagnostics.csv -print
