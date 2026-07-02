#!/usr/bin/env python3
import argparse
import subprocess
import sys
import itertools
import os
import math
import glob
import re

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

def parse_range_or_list(value_str, allow_step=False):
    """
    Parses a string representing a range (e.g. '3-6', with step '10-20:2', '1000-10000:1.5exp' or '1000-10000:2exp')
    or a comma-separated list.
    Returns a list of parsed values (integers or strings).
    """
    if not value_str:
        return []
    
    items = []
    for part in value_str.split(','):
        part = part.strip()
        
        step = 1
        is_exp = False
        multiplier = 2.0
        range_part = part
        if allow_step and ':' in part:
            subparts_step = part.split(':')
            if len(subparts_step) == 2:
                range_part = subparts_step[0].strip()
                step_str = subparts_step[1].strip()
                if step_str.endswith("exp"):
                    is_exp = True
                    mult_str = step_str[:-3].strip()
                    if mult_str:
                        try:
                            multiplier = float(mult_str)
                            if multiplier <= 0.0:
                                print(f"Warning: Invalid multiplier '{multiplier}'. Multiplier must be positive. Using default step 1.", file=sys.stderr)
                                multiplier = 2.0
                        except ValueError:
                            print(f"Warning: Could not parse multiplier '{mult_str}'. Using default multiplier 2.0.", file=sys.stderr)
                            multiplier = 2.0
                    else:
                        multiplier = 2.0
                else:
                    try:
                        step = int(step_str)
                        if step <= 0:
                            print(f"Warning: Invalid step '{step}'. Step must be a positive integer. Using default step 1.", file=sys.stderr)
                            step = 1
                    except ValueError:
                        print(f"Warning: Could not parse step '{step_str}'. Using default step 1.", file=sys.stderr)
                        step = 1

        if '-' in range_part:
            subparts = range_part.split('-')
            if len(subparts) == 2:
                try:
                    start = int(subparts[0].strip())
                    end = int(subparts[1].strip())
                    if is_exp:
                        if start <= 0 or (start > end and end <= 0):
                            print(f"Warning: Exponential/multiplier step is only supported for positive numbers (start={start}, end={end}).", file=sys.stderr)
                            if start <= end:
                                items.extend(range(start, end + 1, 1))
                            else:
                                items.extend(range(start, end - 1, -1))
                        else:
                            curr = start
                            if start <= end:
                                while curr <= end:
                                    items.append(curr)
                                    next_curr = int(round(curr * multiplier))
                                    if next_curr <= curr:
                                        next_curr = curr + 1
                                    curr = next_curr
                            else:
                                while curr >= end:
                                    items.append(curr)
                                    next_curr = int(round(curr * multiplier))
                                    if next_curr >= curr:
                                        next_curr = curr - 1
                                    curr = next_curr
                    else:
                        if start <= end:
                            items.extend(range(start, end + 1, step))
                        else:
                            items.extend(range(start, end - 1, -step))
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

def parse_likwid_output(stdout, group_name):
    """
    Parses the stdout of a likwid-perfctr run.
    Returns a dictionary of {f"{group_name}_{metric/event}": value}.
    """
    results = {}
    lines = stdout.strip().split("\n")
    
    in_raw_table = False
    in_metric_table = False
    raw_header_seen = False
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        # Detect table markers
        if line.startswith("TABLE,"):
            parts = line.split(",")
            if len(parts) >= 4:
                table_type = parts[2]  # e.g., "Group 1 Raw" or "Group 1 Metric"
                if "Raw" in table_type:
                    in_raw_table = True
                    in_metric_table = False
                    raw_header_seen = False
                elif "Metric" in table_type:
                    in_raw_table = False
                    in_metric_table = True
            continue
        
        if in_raw_table:
            if line.startswith("STRUCT,") or line.startswith("TABLE,"):
                in_raw_table = False
                raw_header_seen = False
            else:
                parts = [p.strip() for p in line.split(",")]
                if parts[0] == "Event":
                    raw_header_seen = True
                    continue
                
                if raw_header_seen:
                    if len(parts) >= 3:
                        name = parts[0]
                        val = parts[2]
                        results[f"{group_name}_{name}"] = val
                else:
                    if len(parts) >= 2:
                        name = parts[0]
                        val = parts[1]
                        if name and name != "Region Info":
                            results[f"{group_name}_{name}"] = val
                            
        elif in_metric_table:
            if line.startswith("STRUCT,") or line.startswith("TABLE,"):
                in_metric_table = False
            else:
                parts = [p.strip() for p in line.split(",")]
                if parts[0] == "Metric":
                    continue
                if len(parts) >= 2:
                    name = parts[0]
                    val = parts[1]
                    if name:
                        results[f"{group_name}_{name}"] = val
                        
    return results

def gather_slurm_results(directory, output_file, pretty):
    """Parses all Slurm .out files in the directory and builds the final CSV."""
    header = "ncells,ndofs,fe_degree,threads,mpi_procs,setup_time,assemble_time,solve_time,solver_iteration_time,solver_iter_time_per_dof,iterations,mesh_str"
    results_cache = []

    # Locate all SLURM output files created by this script
    files = glob.glob(os.path.join(directory, "bench_*.out"))

    if not files:
        print(f"No output files found matching 'bench_*.out' in {directory}", file=sys.stderr)
        sys.exit(1)

    for filepath in files:
        # Extract thread and mpi info from filename: bench_c10_d1_mf_t4_p1_12345.out
        filename = os.path.basename(filepath)
        match = re.search(r'bench_c\d+_d\d+_[a-z]+_t(\d+)_p(\d+)_', filename)
        if not match:
            continue

        t_val, mpi_val = match.group(1), match.group(2)

        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split(",")
                # Look for the exact 10-column output generated by the C++ code
                if len(parts) == 10 and parts[0] != "ncells":
                    final_parts = [
                        parts[0], parts[1], parts[2], t_val, mpi_val,
                        parts[4], parts[5], parts[6], parts[7], parts[8],
                        parts[3], parts[9]
                    ]
                    results_cache.append(",".join(final_parts))
                    break

    # Sort results by ncells, then by polynomial degree
    results_cache.sort(key=lambda x: (int(x.split(",")[0]), int(x.split(",")[2])))

    # Output formatting
    out_f = open(output_file, "w") if output_file else sys.stdout
    if pretty:
        headers = header.split(",")
        rows = [row.split(",") for row in results_cache]
        col_widths = [len(h) for h in headers]
        for row in rows:
            for col_idx, val in enumerate(row):
                col_widths[col_idx] = max(col_widths[col_idx], len(val))
        fmt = " | ".join(f"{{:<{w}}}" for w in col_widths)
        print(fmt.format(*headers), file=out_f)
        print("-+-".join("-" * w for w in col_widths), file=out_f)
        for row in rows:
            print(fmt.format(*row), file=out_f)
    else:
        print(header, file=out_f)
        for row in results_cache:
            print(row, file=out_f)

    if output_file:
        out_f.close()
        print(f"Results gathered and saved to {output_file}")

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
    parser.add_argument(
        "--likwid",
        action="store_true",
        help="Use likwid-perfctr to measure performance counters."
    )
    parser.add_argument(
        "--likwid-groups",
        default="MEM_DP",
        help="Comma-separated list of LIKWID performance groups to measure (e.g. MEM_DP,L2,L3). Default: MEM_DP."
    )

    parser.add_argument("--slurm", action="store_true", help="Submit runs as G100 Slurm batch jobs.")
    parser.add_argument("--sif", default="../hpc.sif", help="Path to the Singularity .sif file.")
    parser.add_argument("--account", default="pMI25_DCMIC", help="Slurm account.")
    parser.add_argument("--partition", default="g100_usr_prod", help="Slurm partition.")
    parser.add_argument("--mem", default="150GB", help="Memory allocation for the Slurm job.")
    parser.add_argument("--time", default="01:00:00", help="Max Slurm time limit (HH:MM:SS).")
    parser.add_argument("--gather", metavar="DIR", help="Gather results from Slurm .out files in DIR.")

    args = parser.parse_args()

    # --- MODE 1: GATHER SLURM RESULTS ---
    if args.gather:
        gather_slurm_results(args.gather, args.output, args.pretty)
        sys.exit(0)
    
    degree_values = parse_range_or_list(args.degree)

    if not degree_values:
        print("Error: No polynomial degrees specified or parsed successfully.", file=sys.stderr)
        sys.exit(1)

    # Resolve threads and mpi-procs mutual exclusivity and defaults
    if args.threads and args.mpi_procs:
        print("Error: --threads and --mpi-procs are mutually exclusive.", file=sys.stderr)
        sys.exit(1)

    impl_values = parse_impl(args.impl)

    likwid_groups = []
    if args.likwid_groups:
        likwid_groups = [g.strip() for g in args.likwid_groups.split(",") if g.strip()]

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
            if args.likwid:
                thread_values = [1]
            else:
                thread_values = [os.cpu_count() or 1]

        for t in thread_values:
            if not isinstance(t, int) or t <= 0:
                print(f"Error: Invalid thread count '{t}'. Must be a positive integer.", file=sys.stderr)
                sys.exit(1)

    if args.likwid:
        if any(p > 1 for p in mpi_procs_values):
            print("Error: When using --likwid, you cannot use --mpi-procs.", file=sys.stderr)
            sys.exit(1)
        if any(t > 1 for t in thread_values):
            print("Error: When using --likwid, you cannot use more than 1 thread (-t 1).", file=sys.stderr)
            sys.exit(1)

    binary_path = args.binary
    if not os.path.exists(binary_path):
        print(f"Error: Executable not found at '{binary_path}'. Did you compile the project?", file=sys.stderr)
        sys.exit(1)

    # Determine permutations
    permutations = []
    if args.cells_per_dim:
        c_values = parse_range_or_list(args.cells_per_dim, allow_step=True)
        if not c_values:
            print("Error: No cells-per-dim values specified or parsed successfully.", file=sys.stderr)
            sys.exit(1)
        # Cartesian product including solver implementations and thread counts
        for c_val, deg_val, impl_val, t_val, mpi_val in itertools.product(c_values, degree_values, impl_values, thread_values, mpi_procs_values):
            permutations.append((c_val, deg_val, impl_val, t_val, mpi_val))
    else:
        dofs_values = parse_range_or_list(args.dofs, allow_step=True)
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

    base_cols = ["ncells", "ndofs", "fe_degree", "impl", "threads", "mpi_procs", "setup_time", "assemble_time", "solve_time", "solver_iteration_time", "solver_iter_time_per_dof", "iterations", "mesh_str"]
    base_header = ",".join(base_cols)
    full_cols = base_cols[:]
    header = base_header
    
    header_written = False
    if not args.likwid:
        if out_file:
            out_file.write(base_header + "\n")
            out_file.flush()
        header_written = True

    results_cache = []
    run_groups = likwid_groups if args.likwid else [None]

    for idx, (n_val, deg_val, impl_val, t_val, mpi_val) in enumerate(permutations, 1):
        perm_data = {}
        skip_permutation = False
        
        dofs_est = (n_val * deg_val + 1) ** args.dim if isinstance(n_val, int) else "unknown"

        job_base_name = f"c{n_val}_d{deg_val}_{impl_val}_t{t_val}_p{mpi_val}"
        
        # Build the actual bash command to run the C++ executable
        run_cmd = f"{args.binary} {n_val} {deg_val} {impl_val} {args.iters}"
            
        # SLURM EXECUTION (Cineca Template)
        if args.slurm:
            script_file = f"submit_{job_base_name}.sh"
            
            sing_cmd = f"singularity exec --cleanenv -B $CINECA_SCRATCH {args.sif} bash -c \"source /u/sw/lmod/lmod/init/bash; export MODULEPATH=/u/sw/toolchains/gcc-glibc/11.2.0/modules:/u/sw/modules/toolchains; module load gcc-glibc dealii && export DEAL_II_NUM_THREADS={t_val} && {run_cmd}\""
            
            if mpi_val > 1:
                final_launch = f"srun -n {mpi_val} {sing_cmd}"
            else:
                final_launch = sing_cmd

            lines = [
                "#!/bin/bash",
                f"#SBATCH --job-name=bench_{job_base_name}",
                f"#SBATCH --partition={args.partition}",
                f"#SBATCH --mem={args.mem}",
                f"#SBATCH --error=bench_{job_base_name}_%j.err",
                f"#SBATCH --output=bench_{job_base_name}_%j.out",
                f"#SBATCH --account={args.account}",
                f"#SBATCH --time={args.time}",
                "#SBATCH --nodes=1",
                f"#SBATCH --ntasks-per-node={mpi_val if mpi_val > 1 else 1}",
                f"#SBATCH --cpus-per-task={t_val}",
                "",
                "# 1. Load Host Modules",
                "module purge",
                "module load profile/advanced",
                "module load singularity",
                "",
                "# 2. Temporary Directory Management",
                "export SINGULARITY_TMPDIR=$CINECA_SCRATCH/sing_tmp_$SLURM_JOBID",
                "export SINGULARITY_CACHEDIR=$CINECA_SCRATCH/sing_cache_$SLURM_JOBID",
                "mkdir -p $SINGULARITY_TMPDIR",
                "mkdir -p $SINGULARITY_CACHEDIR",
                "",
                "# 3. Container Execution",
                "# We map /g100_work and $CINECA_SCRATCH, set the threads, and run the binary",
                f"{final_launch}",
                "",
                "# 4. Final Cleanup",
                "rm -rf $SINGULARITY_TMPDIR",
                "rm -rf $SINGULARITY_CACHEDIR"
            ]
            
            with open(script_file, "w") as f:
                f.write("\n".join(lines) + "\n")
                
            print(f"[{idx}/{len(permutations)}] Submitting to Slurm -> {script_file}", file=sys.stderr)
            subprocess.run(["sbatch", script_file], check=True)

        # LOCAL EXECUTION (Original Script Logic)
        else:


            # Display progress info on stderr so it doesn't pollute CSV stdout
            print(f"[{idx}/{len(permutations)}] Running benchmark: cells_per_dim={n_val}, dofs={dofs_est}, degree={deg_val}, solver={impl_val}, threads={t_val}, mpi_procs={mpi_val}, iterations={args.iters}", file=sys.stderr, flush=True)
        
            run_env = os.environ.copy()
            run_env["DEAL_II_NUM_THREADS"] = str(t_val)
        
            for g_idx, group in enumerate(run_groups):
                if args.likwid:
                    cmd = [
                        "likwid-perfctr",
                        "-C", "0",
                        "-g", group,
                        "-m",
                        "-O",
                        binary_path,
                        str(n_val),
                        str(deg_val),
                        impl_val,
                        str(args.iters)
                    ]
                else:
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
            
                try:
                    result = subprocess.run(cmd, env=run_env, capture_output=True, text=True, check=True)
                
                    # Parse stdout to locate the data row
                    lines = result.stdout.strip().split("\n")
                
                    if g_idx == 0:
                        base_parsed = False
                        for line in lines:
                            parts = [p.strip() for p in line.split(",")]
                            if len(parts) == 10 and parts[0] != "ncells":
                                perm_data["ncells"] = parts[0]
                                perm_data["ndofs"] = parts[1]
                                perm_data["fe_degree"] = parts[2]
                                perm_data["impl"] = impl_val
                                perm_data["threads"] = str(t_val)
                                perm_data["mpi_procs"] = str(mpi_val)
                                perm_data["setup_time"] = parts[4]
                                perm_data["assemble_time"] = parts[5]
                                perm_data["solve_time"] = parts[6]
                                perm_data["solver_iteration_time"] = parts[7]
                                perm_data["solver_iter_time_per_dof"] = parts[8]
                                perm_data["iterations"] = parts[3]
                                perm_data["mesh_str"] = parts[9]
                                base_parsed = True
                                break
                        if not base_parsed:
                            print(f"Warning: Could not parse CSV output row from command: {' '.join(cmd)}", file=sys.stderr)
                            if result.stdout:
                                print(f"Output was:\n{result.stdout}", file=sys.stderr)
                
                    if args.likwid:
                        group_data = parse_likwid_output(result.stdout, group)
                        perm_data.update(group_data)
                    
                except subprocess.CalledProcessError as e:
                    print(f"Error running command: {' '.join(cmd)}", file=sys.stderr)
                    print(f"Exit code: {e.returncode}", file=sys.stderr)
                    if e.stderr:
                        print(f"Error output:\n{e.stderr.strip()}", file=sys.stderr)
                    print("-" * 40, file=sys.stderr)
                    skip_permutation = True
                    break
        
            if skip_permutation:
                continue
            
            if args.likwid and not header_written:
                likwid_cols = sorted([k for k in perm_data.keys() if k not in base_cols])
                full_cols = base_cols + likwid_cols
                header = ",".join(full_cols)
                if out_file:
                    out_file.write(header + "\n")
                    out_file.flush()
                header_written = True
            
            row_values = [str(perm_data.get(col, "")) for col in full_cols]
            data_row = ",".join(row_values)
            results_cache.append(data_row)
            if out_file:
                out_file.write(data_row + "\n")
                out_file.flush()

    # Print the clean accumulated CSV results to stdout at the very end
    if args.pretty:
        if not args.slurm:
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
