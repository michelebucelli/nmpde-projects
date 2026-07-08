#!/usr/bin/env python3
import sys
import os
import re

def parse_results(results_file):
    data = []
    pattern = re.compile(r'^(\S+)\s+-\s+(.*?):\s+([0-9.]+)\s+seconds')
    with open(results_file, 'r') as f:
        for line in f:
            match = pattern.match(line.strip())
            if match:
                job_id = match.group(1)
                config = match.group(2).strip()
                # Clean up multiple spaces
                config = re.sub(r'\s+', ' ', config)
                time_s = float(match.group(3))
                data.append((job_id, config, time_s))
    return data

def parse_iterations(out_file):
    iterations = []
    pattern = re.compile(r'Timestep\s+\d+.*:\s*(\d+)\s+\w+\s+iterations')
    if not os.path.isfile(out_file):
        return 0.0
    with open(out_file, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                iterations.append(int(match.group(1)))
    return sum(iterations) / len(iterations) if iterations else 0.0

def text_plot(data):
    if not data:
        print("No data to display.")
        return
        
    # Sort by CPU count
    data.sort(key=lambda x: int(re.match(r'^(\d+)', x[0]).group(1)) if re.match(r'^(\d+)', x[0]) else 0)
        
    print("\n=== Text-Based Performance Plot ===")
    max_len = max(len(d[0]) for d in data)
    
    # Times scale
    max_time = max(d[1] for d in data)
    scale_time = 25.0 / max_time if max_time > 0 else 1.0
    
    # Iters scale
    max_iter = max(d[3] for d in data)
    scale_iter = 25.0 / max_iter if max_iter > 0 else 1.0
    
    for config, mean_time, std_time, mean_iter, std_iter in data:
        bar_time = '#' * int(mean_time * scale_time)
        bar_iter = '*' * int(mean_iter * scale_iter)
        time_str = f"{mean_time:.2f} ± {std_time:.2f}m" if std_time > 0 else f"{mean_time:.2f}m"
        iter_str = f"{mean_iter:.1f} ± {std_iter:.1f}" if std_iter > 0 else f"{mean_iter:.1f}"
        print(f"{config:<{max_len}} | Time: {bar_time:<25} ({time_str}) | Iters: {bar_iter:<25} ({iter_str})")
    print("===================================\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: plot_scaling.py <results_directory>")
        sys.exit(1)
        
    directory = sys.argv[1]
    results_file = os.path.join(directory, "parallel_scalability.txt")
    if not os.path.isfile(results_file):
        print(f"Error: Results file '{results_file}' not found.")
        sys.exit(1)
        
    raw_data = parse_results(results_file)
    if not raw_data:
        print(f"Error: No valid result data found in '{results_file}'.")
        sys.exit(1)
        
    from collections import defaultdict
    import math
    
    # Group times and iterations by configuration
    config_times = defaultdict(list)
    config_iters = defaultdict(list)
    for job_id, config, time_s in raw_data:
        out_file = os.path.join(directory, f"{job_id}.out")
        avg_iter = parse_iterations(out_file)
        config_times[config].append(time_s / 60.0)
        config_iters[config].append(avg_iter)
        
    def get_stats(vals):
        if not vals:
            return 0.0, 0.0
        m = sum(vals) / len(vals)
        if len(vals) > 1:
            var = sum((x - m) ** 2 for x in vals) / (len(vals) - 1)
            sd = math.sqrt(var)
        else:
            sd = 0.0
        return m, sd

    data = []
    cpu_counts = []
    for config in config_times:
        mean_time, std_time = get_stats(config_times[config])
        mean_iter, std_iter = get_stats(config_iters[config])
        data.append((config, mean_time, std_time, mean_iter, std_iter))
        
        # parse cpu count
        match = re.match(r'^(\d+)', config)
        if match:
            cpu_counts.append(int(match.group(1)))
        else:
            cpu_counts.append(1)
            
    # Attempt to plot with matplotlib
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        # Modern styling
        plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
        
        # Sort CPU configurations numerically
        sorted_data = sorted(zip(cpu_counts, data), key=lambda x: x[0])
        sorted_cpus = [x[0] for x in sorted_data]
        sorted_configs = [x[1][0] for x in sorted_data]
        sorted_times = [x[1][1] for x in sorted_data]
        std_times = [x[1][2] for x in sorted_data]
        sorted_iters = [x[1][3] for x in sorted_data]
        std_iters = [x[1][4] for x in sorted_data]
        
        # Setup ticks and labels using only the number of cores in use
        x_ticks = sorted_cpus
        x_labels = [str(c) for c in sorted_cpus]
        
        # Plot 1: Execution Time (Line Plot with Std Dev)
        fig1, ax1 = plt.subplots(figsize=(8, 5))
        ax1.errorbar(sorted_cpus, sorted_times, yerr=std_times, marker='o', markersize=6, linewidth=2, color='#2b5c8f',
                     ecolor='#7f8c8d', capsize=4, elinewidth=1.5, label='Execution Time')
        
        # Add labels above points
        max_time = max(m + s for m, s in zip(sorted_times, std_times)) if sorted_times else 1.0
        ax1.set_ylim(0, max_time * 1.18)
        for x, y, std in zip(sorted_cpus, sorted_times, std_times):
            label_text = f'{y:.2f}±{std:.2f}m' if std > 0 else f'{y:.2f}m'
            ax1.text(x, y + std + (max_time * 0.03), label_text, 
                     ha='center', va='bottom', fontsize=9, fontweight='bold', color='#2c3e50')
        
        ax1.set_xlabel('Number of cores in use', fontsize=12, fontweight='bold', labelpad=10)
        ax1.set_ylabel('Execution Time (minutes)', fontsize=12, fontweight='bold', labelpad=10)
        ax1.set_title('Scaling on parallel infrastructure', fontsize=14, fontweight='bold', pad=15)
        ax1.grid(True, linestyle='--', alpha=0.5)
        ax1.set_xticks(x_ticks)
        ax1.set_xticklabels(x_labels)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_color('#bdc3c7')
        ax1.spines['bottom'].set_color('#bdc3c7')
        plt.tight_layout()
        
        output_image1 = os.path.join(directory, "execution_times.png")
        plt.savefig(output_image1, dpi=300)
        plt.close(fig1)
        print(f"Successfully generated execution times line plot: {output_image1}")
        
        # Plot 2: Average Iterations (Line Plot with Std Dev)
        fig2, ax2 = plt.subplots(figsize=(8, 5))
        ax2.errorbar(sorted_cpus, sorted_iters, yerr=std_iters, marker='s', markersize=6, linewidth=2, color='#e74c3c',
                     ecolor='#7f8c8d', capsize=4, elinewidth=1.5, label='Average Iterations')
        
        # Add labels above points
        max_iter = max(m + s for m, s in zip(sorted_iters, std_iters)) if sorted_iters else 1.0
        ax2.set_ylim(0, max_iter * 1.18)
        for x, y, std in zip(sorted_cpus, sorted_iters, std_iters):
            label_text = f'{y:.1f}±{std:.1f}' if std > 0 else f'{y:.1f}'
            ax2.text(x, y + std + (max_iter * 0.03 if max_iter > 0 else 0.2), label_text, 
                     ha='center', va='bottom', fontsize=9, fontweight='bold', color='#2c3e50')
        
        ax2.set_xlabel('Number of cores in use', fontsize=12, fontweight='bold', labelpad=10)
        ax2.set_ylabel('Average Iterations per Timestep', fontsize=12, fontweight='bold', labelpad=10)
        ax2.set_title('Average Iterations based on MPI proccesses', fontsize=14, fontweight='bold', pad=15)
        ax2.grid(True, linestyle='--', alpha=0.5)
        ax2.set_xticks(x_ticks)
        ax2.set_xticklabels(x_labels)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_color('#bdc3c7')
        ax2.spines['bottom'].set_color('#bdc3c7')
        plt.tight_layout()
        
        output_image2 = os.path.join(directory, "average_iterations.png")
        plt.savefig(output_image2, dpi=300)
        plt.close(fig2)
        print(f"Successfully generated average iterations line plot: {output_image2}")
            
    except ImportError:
        print("Matplotlib is not installed. Generating text-based ASCII plot instead:")
        text_plot(data)
        print("Tip: Install matplotlib (`pip install matplotlib`) to generate high-quality PNG plots.")

if __name__ == '__main__':
    main()
