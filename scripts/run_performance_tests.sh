#!/bin/bash
# Performance monitoring script for Azure deployments
# Collects CPU, memory, disk, and timing metrics

set -e

BUILD_DIR=${1:-build}
REPORT_DIR=${2:-reports}
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
REPORT_FILE="${REPORT_DIR}/performance_${TIMESTAMP}.json"
METRICS_FILE="${REPORT_DIR}/system_metrics_${TIMESTAMP}.txt"

mkdir -p "${REPORT_DIR}"

echo "==============================================="
echo "Performance Test Runner for Azure"
echo "==============================================="
echo "Build directory: ${BUILD_DIR}"
echo "Report directory: ${REPORT_DIR}"
echo "Timestamp: ${TIMESTAMP}"
echo ""

# Collect system information
collect_system_info() {
    cat > "${METRICS_FILE}" << EOF
=== Performance Test Report ===
Date: $(date -Iseconds)
Hostname: $(hostname)
Kernel: $(uname -r)

=== CPU Information ===
$(lscpu 2>/dev/null || echo "lscpu not available")

=== Memory Information ===
$(free -h 2>/dev/null || echo "free command not available")

=== Disk Information ===
$(df -h 2>/dev/null || echo "df command not available")

=== Process Information ===
CPU Count: $(nproc)
Max Processes: $(ulimit -u)
Max Memory: $(ulimit -m)

EOF
}

# Initialize JSON report
init_json_report() {
    cat > "${REPORT_FILE}" << EOF
{
  "timestamp": "$(date -Iseconds)",
  "hostname": "$(hostname)",
  "cpu_count": $(nproc),
  "memory_gb": $(free -g 2>/dev/null | awk '/^Mem:/{print $2}' || echo 0),
  "kernel": "$(uname -r)",
  "tests": []
}
EOF
}

# Run a single test with performance monitoring
run_test() {
    local test_name=$1
    local executable=$2
    local args=$3
    
    echo "Running test: ${test_name}"
    echo "Executable: ${executable}"
    echo "Arguments: ${args}"
    
    # Create a temp file for metrics
    local tmp_output="/tmp/perf_test_${test_name}_$$.txt"
    local time_output="/tmp/time_output_${test_name}_$$.txt"
    
    # Start time
    local start_time=$(date +%s.%N)
    
    # Run with time command if available
    if command -v /usr/bin/time &> /dev/null; then
        /usr/bin/time -v ${executable} ${args} > "${tmp_output}" 2>&1 || true
    else
        ${executable} ${args} > "${tmp_output}" 2>&1 || true
    fi
    
    local exit_code=$?
    local end_time=$(date +%s.%N)
    local duration=$(echo "${end_time} - ${start_time}" | bc)
    
    # Extract metrics from output
    local max_mem=0
    if grep -q "Maximum resident set size" "${tmp_output}"; then
        max_mem=$(grep "Maximum resident set size" "${tmp_output}" | awk '{print $6}')
    fi
    
    local cpu_percent=0
    if grep -q "Percent of CPU" "${tmp_output}"; then
        cpu_percent=$(grep "Percent of CPU" "${tmp_output}" | awk '{print $7}' | sed 's/%//')
    fi
    
    # Append to JSON (manual JSON construction to avoid jq dependency)
    echo "  Test: ${test_name}"
    echo "    Duration: ${duration}s"
    echo "    Max Memory: ${max_mem} KB"
    echo "    CPU: ${cpu_percent}%"
    echo "    Exit Code: ${exit_code}"
    echo ""
    
    # Store results
    cat >> "${REPORT_DIR}/test_${test_name}_${TIMESTAMP}.txt" << EOF
Test: ${test_name}
Duration: ${duration} seconds
Max Memory: ${max_mem} KB
CPU Usage: ${cpu_percent}%
Exit Code: ${exit_code}

=== Test Output ===
$(cat "${tmp_output}")
EOF
    
    # Cleanup
    rm -f "${tmp_output}" "${time_output}"
    
    return ${exit_code}
}

# Run MPI test
run_mpi_test() {
    local processes=$1
    local executable=$2
    
    echo "==============================================="
    echo "Running MPI Test with ${processes} processes"
    echo "==============================================="
    
    if ! command -v mpirun &> /dev/null; then
        echo "WARNING: mpirun not found, skipping MPI test"
        return 1
    fi
    
    local start_time=$(date +%s.%N)
    local tmp_output="/tmp/mpi_test_$$.txt"
    
    mpirun -np ${processes} ${executable} > "${tmp_output}" 2>&1 || true
    local exit_code=$?
    
    local end_time=$(date +%s.%N)
    local duration=$(echo "${end_time} - ${start_time}" | bc)
    
    echo "  MPI Run completed in ${duration}s (exit code: ${exit_code})"
    
    cat >> "${REPORT_DIR}/mpi_test_${processes}proc_${TIMESTAMP}.txt" << EOF
MPI Test with ${processes} processes
Duration: ${duration} seconds
Exit Code: ${exit_code}

=== MPI Output ===
$(cat "${tmp_output}")
EOF
    
    rm -f "${tmp_output}"
    return ${exit_code}
}

# Main execution
main() {
    echo "Collecting system information..."
    collect_system_info
    
    echo "Initializing JSON report..."
    init_json_report
    
    # Check if executable exists
    if [ ! -f "${BUILD_DIR}/space-time-adaptivity" ]; then
        echo "ERROR: Executable not found: ${BUILD_DIR}/space-time-adaptivity"
        exit 1
    fi
    
    # Make sure it's executable
    chmod +x "${BUILD_DIR}/space-time-adaptivity"
    
    # Run single-process test
    echo ""
    echo "==============================================="
    echo "Running Single Process Test"
    echo "==============================================="
    run_test "single_process" "${BUILD_DIR}/space-time-adaptivity" ""
    
    # Run MPI tests with different process counts
    if command -v mpirun &> /dev/null; then
        for np in 2 4 8; do
            if [ ${np} -le $(nproc) ]; then
                run_mpi_test ${np} "${BUILD_DIR}/space-time-adaptivity"
            else
                echo "Skipping MPI test with ${np} processes (only $(nproc) cores available)"
            fi
        done
    fi
    
    echo ""
    echo "==============================================="
    echo "Performance tests completed!"
    echo "==============================================="
    echo "Reports saved to: ${REPORT_DIR}"
    echo "  - JSON Report: ${REPORT_FILE}"
    echo "  - System Metrics: ${METRICS_FILE}"
    echo ""
    
    # Display summary
    if [ -f "${METRICS_FILE}" ]; then
        echo "=== System Metrics Summary ==="
        grep -A 20 "CPU Information" "${METRICS_FILE}" || true
    fi
}

# Run main function
main "$@"
