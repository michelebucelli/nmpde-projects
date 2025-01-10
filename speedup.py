import os
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter, LogLocator

# Initialize dictionaries to store times by process count
matrix_free_times = {}
classic_times = {}

# Base directory
base_dir = 'result'

# Traverse the directories within 'result'
for folder_name in os.listdir(base_dir):
    folder_path = os.path.join(base_dir, folder_name)
    
    # Check if it's a folder and the name is a number (process count)
    if os.path.isdir(folder_path) and folder_name.isdigit():
        process_count = int(folder_name)
        
        # Find the relevant 'stdout_*.txt' file in the folder
        for file_name in os.listdir(folder_path):
            if file_name.startswith("stdout_") and file_name.endswith(".txt"):
                file_path = os.path.join(folder_path, file_name)
                
                # Open and read the file
                with open(file_path, 'r') as file:
                    content = file.read()
                    
                    # Extract matrix-free and classic times using regex
                    matrix_free_match = re.search(r"matrix-free simulation with \d+ processors.*?Time:\s+(\d+)\s+ms", content, re.DOTALL)
                    classic_match = re.search(r"classic simulation with \d+ processors.*?Time:\s+(\d+)\s+ms", content, re.DOTALL)
                    
                    if matrix_free_match and classic_match:
                        matrix_free_time = int(matrix_free_match.group(1))
                        classic_time = int(classic_match.group(1))
                        
                        # Store times by process count
                        matrix_free_times[process_count] = matrix_free_time
                        classic_times[process_count] = classic_time

# Define baseline (single-core times) for speedup calculation
baseline_matrix_free = matrix_free_times.get(1)
baseline_classic = classic_times.get(1)

# Calculate speedup
matrix_free_speedup = {p: baseline_matrix_free / t for p, t in matrix_free_times.items()}
classic_speedup = {p: baseline_classic / t for p, t in classic_times.items()}

# Sort the process counts for plotting
sorted_process_counts = sorted(matrix_free_speedup.keys())

# Create a figure with 1 row and 2 columns of subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Plot speedup on the left with logarithmic scale
ax1.plot(sorted_process_counts, [matrix_free_speedup[p] for p in sorted_process_counts], 
         marker='o', label='Matrix-Free Simulation Speedup', color='b')
ax1.plot(sorted_process_counts, [classic_speedup[p] for p in sorted_process_counts], 
         marker='s', label='Classic Simulation Speedup', color='r')
ax1.plot(sorted_process_counts, sorted_process_counts, '--', label='Ideal Speedup', color='gray')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('Number of Processors')
ax1.set_ylabel('Speedup')
ax1.set_title('Speedup of Matrix-Free vs Classic Simulation')
ax1.legend()
ax1.grid(True, 'both', linestyle='-')

# Plot simulation times on the right with logarithmic scale
ax2.plot(sorted_process_counts, [matrix_free_times[p] for p in sorted_process_counts], 
         marker='o', label='Matrix-Free Simulation Time', color='b')
ax2.plot(sorted_process_counts, [classic_times[p] for p in sorted_process_counts], 
         marker='s', label='Classic Simulation Time', color='r')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Number of Processors')
ax2.set_ylabel('Simulation Time (ms)')
ax2.set_title('Simulation Time of Matrix-Free vs Classic Simulation')
ax2.legend()
ax2.grid(True, 'both', linestyle='-')

# Use ScalarFormatter and disable scientific notation for ax1 and ax2 x-axis
formatter = LogFormatter(2)
ax1.xaxis.set_major_formatter(formatter)
ax2.xaxis.set_major_formatter(formatter)
ax1.yaxis.set_major_formatter(formatter)
ax1.yaxis.set_major_formatter(formatter)
ax1.xaxis.set_minor_formatter(formatter)
ax2.xaxis.set_minor_formatter(formatter)

# Set x-axis ticks to only the processor counts where you have data
locator = LogLocator(2)
ax1.xaxis.set_major_locator(locator)
ax1.xaxis.set_minor_locator(locator)
ax2.xaxis.set_major_locator(locator)
ax2.xaxis.set_minor_locator(locator)

# Adjust layout and show plot
plt.tight_layout()
plt.show()
