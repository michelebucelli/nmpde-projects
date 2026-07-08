#!/usr/bin/env pvpython
"""
Script to extract activation times along the line (0,0,0) -> (20,7,3)
for time integration methods (explicit_euler, crank_nicolson, implicit_euler)
and generate a comparison plot.
"""

import sys
import os
import glob
import re
import argparse
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def find_latest_euler_dir(base_dir="build"):
    # Look for both *_euler_comparison and *_euler
    candidates = glob.glob(os.path.join(base_dir, "*_euler_comparison")) + glob.glob(os.path.join(base_dir, "*_euler"))
    if not candidates:
        candidates = glob.glob("*_euler_comparison") + glob.glob("*_euler")
    if not candidates:
        return None
    return max(candidates, key=os.path.getmtime)

def extract_activation_time_over_line(pvtu_file):
    from vtkmodules.vtkIOXML import vtkXMLPUnstructuredGridReader, vtkXMLUnstructuredGridReader
    
    if pvtu_file.endswith('.pvtu'):
        reader = vtkXMLPUnstructuredGridReader()
    else:
        reader = vtkXMLUnstructuredGridReader()
        
    reader.SetFileName(pvtu_file)
    reader.Update()
    grid = reader.GetOutput()
    
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkPolyData
    from vtkmodules.vtkFiltersCore import vtkProbeFilter
    
    # Generate 200 points along the line (0,0,0) -> (20,7,3)
    num_points = 200
    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([20.0, 7.0, 3.0])
    line_pts = np.linspace(p1, p2, num_points)
    
    vtk_pts = vtkPoints()
    for p in line_pts:
        vtk_pts.InsertNextPoint(p)
        
    poly_data = vtkPolyData()
    poly_data.SetPoints(vtk_pts)
    
    probe = vtkProbeFilter()
    probe.SetInputData(poly_data)
    probe.SetSourceData(grid)
    probe.Update()
    
    probed_data = probe.GetOutput()
    arr = probed_data.GetPointData().GetArray('Activation Time')
    if not arr:
        return None, None
        
    from vtkmodules.util.numpy_support import vtk_to_numpy
    act_time = vtk_to_numpy(arr)
    
    # Calculate distance along the line
    pts_out = vtk_to_numpy(probed_data.GetPoints().GetData())
    dist = np.linalg.norm(pts_out, axis=1)
    
    return dist, act_time

def main():
    parser = argparse.ArgumentParser(description="Plot activation time over line (0,0,0) -> (20,7,3) for Euler comparison.")
    parser.add_argument("euler_dir", nargs="?", default=None, help="Path to the euler comparison parent directory")
    args = parser.parse_args()
    
    euler_dir = args.euler_dir
    if not euler_dir:
        euler_dir = find_latest_euler_dir()
        
    if not euler_dir or not os.path.exists(euler_dir):
        print(f"Error: Could not locate euler comparison directory: {euler_dir}")
        sys.exit(1)
        
    print(f"Processing euler comparison directory: {euler_dir}")
    
    methods = ["explicit_euler", "crank_nicolson", "implicit_euler"]
    method_dirs = {}
    
    # Locate subdirectories inside the parent euler_dir
    subdirs = glob.glob(os.path.join(euler_dir, "*"))
    for sd in subdirs:
        if not os.path.isdir(sd):
            continue
        basename = os.path.basename(sd)
        for m in methods:
            if basename.startswith(m):
                method_dirs[m] = sd
                
    if not method_dirs:
        print("Error: No method subdirectories found to process.")
        sys.exit(1)
        
    print("Processing method subdirectories:")
    for m, d in method_dirs.items():
        print(f"  {m:<15} -> {d}")
        
    all_line_data = {}
    for m, d in method_dirs.items():
        pvtu_files = glob.glob(os.path.join(d, "*.pvtu"))
        if not pvtu_files:
            all_vtu = glob.glob(os.path.join(d, "*.vtu"))
            pvtu_files = [f for f in all_vtu if not re.search(r'\.[0-9]{2,}\.vtu$', f)]
            
        if not pvtu_files:
            print(f"Warning: No valid VTU/PVTU files found in {d}")
            continue
            
        try:
            dist, act_time = extract_activation_time_over_line(pvtu_files[0])
            if dist is not None and act_time is not None:
                all_line_data[m] = (dist, act_time)
        except Exception as e:
            print(f"Error processing {m}: {e}")
            
    if not all_line_data:
        print("Error: No line data successfully extracted.")
        sys.exit(1)
        
    # Plot results
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, ax = plt.subplots(figsize=(10, 6))
    
    method_colors = {
        'explicit_euler': '#e74c3c',   # Red
        'crank_nicolson': '#2b5c8f',   # Blue
        'implicit_euler': '#2ecc71'    # Green
    }
    
    method_labels = {
        'explicit_euler': 'Explicit Euler (theta = 0.0)',
        'crank_nicolson': 'Crank-Nicolson (theta = 0.5)',
        'implicit_euler': 'Implicit Euler (theta = 1.0)'
    }
    
    for m, (dist, act_time) in sorted(all_line_data.items()):
        color = method_colors.get(m, '#7f8c8d')
        label = method_labels.get(m, m)
        
        valid = ~np.isnan(act_time)
        if np.any(valid):
            ax.plot(dist[valid], act_time[valid], label=label, color=color, linewidth=2.5)
            
    ax.set_xlabel('Distance along line (0,0,0) -> (20,7,3) [mm]', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_ylabel('Activation Time [ms]', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_title('Activation Time Profile comparison for Time Discretization Methods\nLine (0,0,0) to (20,7,3)', fontsize=13, fontweight='bold', pad=15)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(fontsize=10, frameon=True, facecolor='white', edgecolor='#bdc3c7')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#bdc3c7')
    ax.spines['bottom'].set_color('#bdc3c7')
    
    plt.tight_layout()
    output_png = os.path.join(euler_dir, "euler_comparison_over_line.png")
    plt.savefig(output_png, dpi=300)
    plt.close(fig)
    print(f"\n[Success] Comparison plot saved successfully to: {output_png}")

    # Plot 2: Zoomed Plot (x >= 20.5)
    fig_zoom, ax_zoom = plt.subplots(figsize=(10, 6))
    for m, (dist, act_time) in sorted(all_line_data.items()):
        color = method_colors.get(m, '#7f8c8d')
        label = method_labels.get(m, m)
        
        valid = ~np.isnan(act_time)
        if np.any(valid):
            ax_zoom.plot(dist[valid], act_time[valid], label=label, color=color, linewidth=2.5)
            
    ax_zoom.set_xlabel('Distance along line (0,0,0) -> (20,7,3) [mm]', fontsize=12, fontweight='bold', labelpad=10)
    ax_zoom.set_ylabel('Activation Time [ms]', fontsize=12, fontweight='bold', labelpad=10)
    ax_zoom.set_title('Activation Time Profile (Zoomed x >= 20.5 mm)\nLine (0,0,0) to (20,7,3)', fontsize=13, fontweight='bold', pad=15)
    ax_zoom.grid(True, linestyle='--', alpha=0.5)
    ax_zoom.legend(fontsize=10, frameon=True, facecolor='white', edgecolor='#bdc3c7')
    
    # Calculate global max distance
    max_dist = 21.4
    for dist, _ in all_line_data.values():
        if len(dist) > 0:
            max_dist = max(max_dist, np.max(dist))
            
    ax_zoom.set_xlim(20.5, max_dist)
    
    # Find min/max y-values in the zoomed range to set nice y-limits
    y_vals_in_zoom = []
    for dist, act_time in all_line_data.values():
        valid_idx = (dist >= 20.5) & (~np.isnan(act_time))
        if np.any(valid_idx):
            y_vals_in_zoom.extend(act_time[valid_idx])
            
    if y_vals_in_zoom:
        ymin, ymax = min(y_vals_in_zoom), max(y_vals_in_zoom)
        padding = (ymax - ymin) * 0.1 if ymax > ymin else 1.0
        ax_zoom.set_ylim(ymin - padding, ymax + padding)
        
    ax_zoom.spines['top'].set_visible(False)
    ax_zoom.spines['right'].set_visible(False)
    ax_zoom.spines['left'].set_color('#bdc3c7')
    ax_zoom.spines['bottom'].set_color('#bdc3c7')
    
    plt.tight_layout()
    output_png_zoom = os.path.join(euler_dir, "euler_comparison_over_line_zoom.png")
    plt.savefig(output_png_zoom, dpi=300)
    plt.close(fig_zoom)
    print(f"[Success] Zoomed comparison plot saved successfully to: {output_png_zoom}")

    # Plot 3: Combined Main Plot with Zoomed Inset
    fig_inset, ax_main = plt.subplots(figsize=(10, 6))
    
    # 1. Plot main curves on ax_main
    for m, (dist, act_time) in sorted(all_line_data.items()):
        color = method_colors.get(m, '#7f8c8d')
        label = method_labels.get(m, m)
        valid = ~np.isnan(act_time)
        if np.any(valid):
            ax_main.plot(dist[valid], act_time[valid], label=label, color=color, linewidth=2.5)
            
    ax_main.set_xlabel('Distance along line (0,0,0) -> (20,7,3) [mm]', fontsize=12, fontweight='bold', labelpad=10)
    ax_main.set_ylabel('Activation Time [ms]', fontsize=12, fontweight='bold', labelpad=10)
    ax_main.set_title('Activation Time Profile with Zoomed Inset\nLine (0,0,0) to (20,7,3)', fontsize=13, fontweight='bold', pad=15)
    ax_main.grid(True, linestyle='--', alpha=0.5)
    ax_main.legend(fontsize=10, frameon=True, facecolor='white', edgecolor='#bdc3c7')
    
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)
    ax_main.spines['left'].set_color('#bdc3c7')
    ax_main.spines['bottom'].set_color('#bdc3c7')
    
    # 2. Add inset axes
    # Position: [x, y, width, height] as fractions of the main axes size
    ax_ins = ax_main.inset_axes([0.55, 0.15, 0.35, 0.35])
    
    # Plot curves on ax_ins
    for m, (dist, act_time) in sorted(all_line_data.items()):
        color = method_colors.get(m, '#7f8c8d')
        valid = ~np.isnan(act_time)
        if np.any(valid):
            ax_ins.plot(dist[valid], act_time[valid], color=color, linewidth=2.0)
            
    # Set inset limits
    max_dist = 21.4
    for dist, _ in all_line_data.values():
        if len(dist) > 0:
            max_dist = max(max_dist, np.max(dist))
            
    ax_ins.set_xlim(20.5, max_dist)
    
    # Set inset y-limits based on zoomed values
    y_vals_in_zoom = []
    for dist, act_time in all_line_data.values():
        valid_idx = (dist >= 20.5) & (~np.isnan(act_time))
        if np.any(valid_idx):
            y_vals_in_zoom.extend(act_time[valid_idx])
            
    if y_vals_in_zoom:
        ymin, ymax = min(y_vals_in_zoom), max(y_vals_in_zoom)
        padding = (ymax - ymin) * 0.1 if ymax > ymin else 1.0
        ax_ins.set_ylim(ymin - padding, ymax + padding)
        
    ax_ins.grid(True, linestyle=':', alpha=0.6)
    
    # Draw connections indicating the zoom region
    try:
        ax_main.indicate_inset_zoom(ax_ins, edgecolor="black", alpha=0.3)
    except Exception:
        pass
        
    plt.tight_layout()
    output_png_combined = os.path.join(euler_dir, "euler_comparison_combined.png")
    plt.savefig(output_png_combined, dpi=300)
    plt.close(fig_inset)
    print(f"[Success] Combined comparison plot with inset saved successfully to: {output_png_combined}")

if __name__ == "__main__":
    main()
