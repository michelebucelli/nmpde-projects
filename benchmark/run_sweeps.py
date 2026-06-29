#!/usr/bin/env python3
import argparse
import subprocess
import sys
import itertools
import os
import math

def calculate_n_from_dofs(min_dofs, degree, dim):
    """
    Calculates the minimum cells per dimension N required to have at least min_dofs.
    Formula: (N * degree + 1)^dim >= min_dofs
    """
    if min_dofs <= 0:
        return 1
    val = (min_dofs ** (1.0 / dim)) - 1.0
    N = int(math.ceil(val / degree))
    return max(1, N)

def parse_range_or_list(value_str):
    """
    Parses a string representing a range (e.g. '3-6') or a comma-separated list (e.g. '3,5,7').
    Returns a list of parsed values (integers or strings).
    """
    if not value_str:
        return []
    
    items = []
    for part in value_str.split(','):
        part = part.strip()
        if '-' in part:
            subparts = part.split('-')
            if len(subparts) == 2:
                try:
                    start = int(subparts[0])
                    end = int(subparts[1])
                    if start <= end:
                        items.extend(range(start, end + 1))
                    else:
                        items.extend(range(start, end - 1, -1))
                    continue
                except ValueError:
                    pass
        try:
            items.append(int(part))
        except ValueError:
            items.append(part)
    return items

def parse_impl(impl_str):
    """
    Parses comma-separated solver implementations. Supports 'mb', 'mf', or 'both'.
    """
    if not impl_str:
        return ["mf"]
    impls = [x.strip() for x in impl_str.split(",")]
    resolved = []
    for impl in impls:
        if impl == "both":
            resolved.extend(["mb", "mf"])
        elif impl in ["mb", "mf"]:
            resolved.append(impl)
        else:
            print(f"Error: Invalid solver implementation '{impl}'. Must be 'mb', 'mf', or 'both'.", file=sys.stderr)
            sys.exit(1)
    seen = set()
    return [x for x in resolved if not (x in seen or seen.add(x))]

def main():
    parser = argparse.ArgumentParser(
        description="Run sweeps for the single-node Poisson/ADR solver benchmark."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-c", "--cells-per-dim",
        help="Subdivisions (cells per dimension) range/list (e.g., '3-6', '3,5,7') or mesh file path."
    )
    group.add_argument(
        "--dofs",
        help="Target minimum degrees of freedom range/list (e.g., '1000-5000', '1000,5000'). Cells per dimension will be calculated per degree."
    )
    parser.add_argument(
        "-d", "--degree",
        required=True,
        help="Polynomial degree range/list (e.g., '1-4', '1,2,3')."
    )
    parser.add_argument(
        "--dim",
        type=int,
        default=2,
        help="Dimension of the problem (e.g., 2 or 3). Default: 2."
    )
    parser.add_argument(
        "-i", "--impl",
        default="mf",
        help="Solver implementation range/list: 'mf', 'mb', or 'both'. Default: 'mf'."
    )
    parser.add_argument(
        "--iters",
        type=int,
        default=100,
        help="Number of iterations to run the solver for. Default: 100."
    )
    parser.add_argument(
        "--binary",
        default="./build/single-node",
        help="Path to the single-node executable. Default: './build/single-node'."
    )
    parser.add_argument(
        "-t", "--threads",
        help="Number of threads to use (sets DEAL_II_NUM_THREADS). Can be a range/list (e.g. '1-4', '1,2,4'). Mutually exclusive with --mpi-procs."
    )
    parser.add_argument(
        "-p", "--mpi-procs",
        help="Number of MPI processes to use (runs via mpirun). Can be a range/list (e.g. '1-4', '1,2,4'). Default: 1. Mutually exclusive with --threads."
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Print final terminal stdout results in a formatted table instead of CSV."
    )
    parser.add_argument(
        "-o", "--output",
        help="Path to save the output CSV file. If not specified, results are printed to stdout."
    )

    args = parser.parse_args()

    degree_values = parse_range_or_list(args.degree)

    if not degree_values:
        print("Error: No polynomial degrees specified or parsed successfully.", file=sys.stderr)
        sys.exit(1)

    # Resolve threads and mpi-procs mutual exclusivity and defaults
    if args.threads and args.mpi_procs:
        print("Error: --threads and --mpi-procs are mutually exclusive.", file=sys.stderr)
        sys.exit(1)

    impl_values = parse_impl(args.impl)

    if args.mpi_procs:
        mpi_procs_values = parse_range_or_list(args.mpi_procs)
        for p in mpi_procs_values:
            if not isinstance(p, int) or p <= 0:
                print(f"Error: Invalid MPI process count '{p}'. Must be a positive integer.", file=sys.stderr)
                sys.exit(1)
        thread_values = [1]
    else:
        mpi_procs_values = [1]
        if args.threads:
            thread_values = parse_range_or_list(args.threads)
        else:
            thread_values = [os.cpu_count() or 1]

        for t in thread_values:
            if not isinstance(t, int) or t <= 0:
                print(f"Error: Invalid thread count '{t}'. Must be a positive integer.", file=sys.stderr)
                sys.exit(1)

    binary_path = args.binary
    if not os.path.exists(binary_path):
        print(f"Error: Executable not found at '{binary_path}'. Did you compile the project?", file=sys.stderr)
        sys.exit(1)

    # Determine permutations
    permutations = []
    if args.cells_per_dim:
        c_values = parse_range_or_list(args.cells_per_dim)
        if not c_values:
            print("Error: No cells-per-dim values specified or parsed successfully.", file=sys.stderr)
            sys.exit(1)
        # Cartesian product including solver implementations and thread counts
        for c_val, deg_val, impl_val, t_val, mpi_val in itertools.product(c_values, degree_values, impl_values, thread_values, mpi_procs_values):
            permutations.append((c_val, deg_val, impl_val, t_val, mpi_val))
    else:
        dofs_values = parse_range_or_list(args.dofs)
        if not dofs_values:
            print("Error: No DoFs values specified or parsed successfully.", file=sys.stderr)
            sys.exit(1)
        for dofs in dofs_values:
            for deg in degree_values:
                for impl_val in impl_values:
                    for t_val in thread_values:
                        for mpi_val in mpi_procs_values:
                            if isinstance(dofs, int):
                                n_val = calculate_n_from_dofs(dofs, deg, args.dim)
                            else:
                                n_val = dofs
                            permutations.append((n_val, deg, impl_val, t_val, mpi_val))

    # Prepare output stream
    out_file = None
    if args.output:
        try:
            out_file = open(args.output, "w", encoding="utf-8")
        except IOError as e:
            print(f"Error opening output file: {e}", file=sys.stderr)
            sys.exit(1)

    # Print CSV Header
    header = "ncells,ndofs,fe_degree,threads,mpi_procs,setup_time,assemble_time,solve_time,solver_iteration_time,solver_iter_time_per_dof,iterations,mesh_str"
    if out_file:
        out_file.write(header + "\n")
        out_file.flush()

    results_cache = []
    for idx, (n_val, deg_val, impl_val, t_val, mpi_val) in enumerate(permutations, 1):
        if mpi_val > 1:
            cmd = [
                "mpirun",
                "-n", str(mpi_val),
                binary_path,
                str(n_val),
                str(deg_val),
                impl_val,
                str(args.iters)
            ]
        else:
            cmd = [
                binary_path,
                str(n_val),
                str(deg_val),
                impl_val,
                str(args.iters)
            ]
        
        # Display progress info on stderr so it doesn't pollute CSV stdout
        print(f"[{idx}/{len(permutations)}] Running benchmark: cells_per_dim={n_val}, degree={deg_val}, solver={impl_val}, threads={t_val}, mpi_procs={mpi_val}, iterations={args.iters}", file=sys.stderr, flush=True)
        
        # Configure environment variables (e.g. DEAL_II_NUM_THREADS)
        run_env = os.environ.copy()
        run_env["DEAL_II_NUM_THREADS"] = str(t_val)
 
        try:
            result = subprocess.run(cmd, env=run_env, capture_output=True, text=True, check=True)
            
            # Parse result output to locate the data row
            lines = result.stdout.strip().split("\n")
            data_row = None
            for line in lines:
                parts = line.split(",")
                if len(parts) == 10 and parts[0] != "ncells":
                    final_parts = [
                        parts[0],       # ncells
                        parts[1],       # ndofs
                        parts[2],       # fe_degree
                        str(t_val),     # threads
                        str(mpi_val),   # mpi_procs
                        parts[4],       # setup_time
                        parts[5],       # assemble_time
                        parts[6],       # solve_time
                        parts[7],       # solver_iteration_time
                        parts[8],       # solver_iter_time_per_dof
                        parts[3],       # iterations
                        parts[9]        # mesh_str
                    ]
                    data_row = ",".join(final_parts)
                    break
            
            if data_row:
                results_cache.append(data_row)
                if out_file:
                    out_file.write(data_row + "\n")
                    out_file.flush()
            else:
                print(f"Warning: Could not parse CSV output row from command: {' '.join(cmd)}", file=sys.stderr)
                if result.stdout:
                    print(f"Output was:\n{result.stdout}", file=sys.stderr)
                    
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {' '.join(cmd)}", file=sys.stderr)
            print(f"Exit code: {e.returncode}", file=sys.stderr)
            if e.stderr:
                print(f"Error output:\n{e.stderr.strip()}", file=sys.stderr)
            print("-" * 40, file=sys.stderr)

    # Print the clean accumulated CSV results to stdout at the very end
    if args.pretty:
        headers = header.split(",")
        rows = [row.split(",") for row in results_cache]
        
        # Calculate column widths
        col_widths = [len(h) for h in headers]
        for row in rows:
            for col_idx, val in enumerate(row):
                if col_idx < len(col_widths):
                    col_widths[col_idx] = max(col_widths[col_idx], len(val))
        
        fmt = " | ".join(f"{{:<{w}}}" for w in col_widths)
        print(fmt.format(*headers), flush=True)
        print("-+-".join("-" * w for w in col_widths), flush=True)
        for row in rows:
            print(fmt.format(*row), flush=True)
    else:
        print(header, flush=True)
        for row in results_cache:
            print(row, flush=True)

    if out_file:
        out_file.close()

if __name__ == "__main__":
    main()
