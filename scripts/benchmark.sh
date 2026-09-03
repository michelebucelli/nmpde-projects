#!/usr/bin/env bash
# Runs the three benchmark suites of the project and writes one CSV per suite
# into results/. Every row is one solver run; the columns are taken from the
# SUMMARY and VMULT lines the solver prints.
#
#   ./scripts/benchmark.sh [results directory]

set -u

readonly ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly SOLVER="${ROOT}/build/hybrid_matrix_free_solver"
readonly INPUT="${ROOT}/input.prm"
readonly OUT="${1:-${ROOT}/results}"

readonly HEADER="backend,precond,dim,degree,refine,ranks,threads,dofs,iters,setup,rhs,solve,memory_mb,vmult_time,mdofs_per_s,l2,h1"

if [[ ! -x "${SOLVER}" ]]; then
  echo "no solver at ${SOLVER}; build it first" >&2
  exit 1
fi

mkdir -p "${OUT}"

# Pulls "key=value" tokens out of the SUMMARY and VMULT lines and prints them
# in the fixed column order above. A missing key becomes an empty field.
to_csv() {
  awk '
    /^SUMMARY/ || /^VMULT/ {
      for (i = 2; i <= NF; ++i) {
        split($i, kv, "=")
        v[kv[1]] = kv[2]
      }
    }
    END {
      if (v["dofs"] == "") exit 1
      printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
        v["backend"], v["precond"], v["dim"], v["degree"], v["refine"],
        v["ranks"], v["threads"], v["dofs"], v["iters"], v["setup"],
        v["rhs"], v["solve"], v["memory_mb"], v["time"], v["mdofs_per_s"],
        v["l2"], v["h1"]
    }'
}

# run <csv> <ranks> <solver arguments...>
run() {
  local csv="$1" ranks="$2"
  shift 2

  local output
  if ! output="$(mpirun -n "${ranks}" "${SOLVER}" -i "${INPUT}" \
      --set Output/Verbose=false "$@" 2>&1)"; then
    echo "  FAILED: ranks=${ranks} $*" >&2
    echo "${output}" | tail -3 >&2
    return
  fi

  local row
  row="$(echo "${output}" | to_csv)" || { echo "  no result: $*" >&2; return; }
  echo "${row}" >> "${csv}"
  echo "  ${row}"
}

# 1. Matrix-based against matrix-free, at equal problem size. The interesting
#    axis is the polynomial degree and the dimension, since matrix-free wins
#    where the assembled matrix has the most entries per row.
suite_backends() {
  local csv="${OUT}/backends.csv"
  echo "${HEADER}" > "${csv}"
  echo "== backends"

  for backend in matrix_free matrix_based; do
    for spec in "2 1 10" "2 2 9" "2 3 8" "3 1 6" "3 2 5" "3 3 4"; do
      read -r dim degree refine <<< "${spec}"
      run "${csv}" 1 \
        --set "Solver/Backend=${backend}" \
        --set "Problem/Dimension=${dim}" \
        --set "Problem/Degree=${degree}" \
        --refine "${refine}" \
        --threads 1 \
        --set "Output/Benchmark applications=50"
    done
  done
}

# 2. Pure MPI strong scaling of the matrix-free backend: one thread per rank,
#    fixed global problem size.
suite_mpi() {
  local csv="${OUT}/mpi_scaling.csv"
  echo "${HEADER}" > "${csv}"
  echo "== mpi strong scaling"

  for ranks in 1 2 4 8; do
    run "${csv}" "${ranks}" \
      --refine 6 --set Problem/Dimension=3 --set Problem/Degree=2 \
      --threads 1 --set "Output/Benchmark applications=20"
  done
}

# 3. The hybrid grid: every combination of ranks and threads per rank, so the
#    MPI only and thread only rows are the edges of the same table.
suite_hybrid() {
  local csv="${OUT}/hybrid_scaling.csv"
  echo "${HEADER}" > "${csv}"
  echo "== hybrid scaling"

  for ranks in 1 2 4 8; do
    for threads in 1 2 4 8; do
      if (( ranks * threads > 8 )); then
        continue # never oversubscribe the machine
      fi
      run "${csv}" "${ranks}" \
        --refine 6 --set Problem/Dimension=3 --set Problem/Degree=2 \
        --threads "${threads}" --set "Output/Benchmark applications=20"
    done
  done
}

# 4. Convergence of both backends on the same meshes. Equal errors are the
#    cross check that the two implementations discretise the same problem.
suite_convergence() {
  local csv="${OUT}/convergence.csv"
  echo "${HEADER}" > "${csv}"
  echo "== convergence"

  for backend in matrix_free matrix_based; do
    for refine in 3 4 5 6 7; do
      run "${csv}" 1 \
        --set "Solver/Backend=${backend}" --set Problem/Degree=2 \
        --refine "${refine}" --threads 1
    done
  done
}

case "${2:-all}" in
  convergence) suite_convergence ;;
  backends) suite_backends ;;
  mpi)      suite_mpi ;;
  hybrid)   suite_hybrid ;;
  *)        suite_convergence; suite_backends; suite_mpi; suite_hybrid ;;
esac

echo
echo "wrote ${OUT}/*.csv"
