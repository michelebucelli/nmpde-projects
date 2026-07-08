#!/usr/bin/env pvpython
"""
Script to extract activation times from convergence sweep simulations
and generate a 3D plot with Delta t on X-axis, h on Y-axis, and Activation Time on Z-axis.
"""

import sys
import os
os.environ['DISPLAY'] = '' # Force headless offscreen rendering to prevent window flashing
import glob
import re
import argparse
import numpy as np

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for headless environments
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D

# No paraview dependency needed, pure VTK reader used below

def find_latest_convergence_dir(base_dir="build"):
    candidates = glob.glob(os.path.join(base_dir, "*_convergence"))
    if not candidates:
        candidates = glob.glob("*_convergence")
    if not candidates:
        return None
    return max(candidates, key=os.path.getmtime)

def extract_h_dt_from_folder(folder_path):
    folder_name = os.path.basename(folder_path)
    m = re.search(r'h([0-9]+(?:_[0-9]+)?).*?dt([0-9]+(?:_[0-9]+)?)', folder_name)
    if m:
        h_val = float(m.group(1).replace('_', '.'))
        dt_val = float(m.group(2).replace('_', '.'))
        return h_val, dt_val
        
    out_files = glob.glob(os.path.join(folder_path, "*.out"))
    if out_files:
        with open(out_files[0], 'r') as f:
            content = f.read()
            m_h = re.search(r'Mesh size:\s+([0-9.]+)', content)
            m_dt = re.search(r'Delta t:\s+([0-9.]+)', content)
            if m_h and m_dt:
                return float(m_h.group(1)), float(m_h.group(1))
    return None, None

def read_max_activation_time(pvtu_files):
    from vtkmodules.vtkIOXML import vtkXMLPUnstructuredGridReader, vtkXMLUnstructuredGridReader
    
    if isinstance(pvtu_files, list):
        pvtu_file = pvtu_files[0]
    else:
        pvtu_file = pvtu_files
        
    if pvtu_file.endswith('.pvtu'):
        reader = vtkXMLPUnstructuredGridReader()
    else:
        reader = vtkXMLUnstructuredGridReader()
        
    reader.SetFileName(pvtu_file)
    reader.Update()
    grid = reader.GetOutput()
    arr = grid.GetPointData().GetArray('Activation Time')
    if arr:
        # GetRange(0) returns tuple (min, max) for component 0
        return arr.GetRange(0)[1]
    return 0.0

def interpolate_2d_grid(unique_x, unique_y, Z, num_x=100, num_y=100):
    x_fine = np.linspace(unique_x[0], unique_x[-1], num_x)
    y_fine = np.linspace(unique_y[0], unique_y[-1], num_y)
    X_fine, Y_fine = np.meshgrid(x_fine, y_fine, indexing='ij')
    Z_fine = np.zeros_like(X_fine)
    
    for i in range(num_x):
        for j in range(num_y):
            xf = x_fine[i]
            yf = y_fine[j]
            ix = max(0, min(np.searchsorted(unique_x, xf) - 1, len(unique_x) - 2))
            iy = max(0, min(np.searchsorted(unique_y, yf) - 1, len(unique_y) - 2))
            
            x0, x1 = unique_x[ix], unique_x[ix+1]
            y0, y1 = unique_y[iy], unique_y[iy+1]
            
            tx = (xf - x0) / (x1 - x0) if x1 > x0 else 0.0
            ty = (yf - y0) / (y1 - y0) if y1 > y0 else 0.0
            
            z00 = Z[ix, iy]
            z10 = Z[ix+1, iy]
            z01 = Z[ix, iy+1]
            z11 = Z[ix+1, iy+1]
            
            Z_fine[i, j] = (1 - tx) * (1 - ty) * z00 + tx * (1 - ty) * z10 + (1 - tx) * ty * z01 + tx * ty * z11
            
    return x_fine, y_fine, X_fine, Y_fine, Z_fine

def main():
    parser = argparse.ArgumentParser(description="Generate 3D plot of Activation Time vs Delta t and Mesh Size h.")
    parser.add_argument("sweep_dir", nargs="?", default=None, help="Path to the convergence sweep directory")
    args = parser.parse_args()

    sweep_dir = args.sweep_dir
    if not sweep_dir:
        sweep_dir = find_latest_convergence_dir()
        
    if not sweep_dir or not os.path.exists(sweep_dir):
        print(f"Error: Could not locate convergence sweep directory: {sweep_dir}")
        sys.exit(1)

    print(f"Processing convergence sweep directory: {sweep_dir}")

    subdirs = sorted(glob.glob(os.path.join(sweep_dir, "*")))
    data_dict = {}

    for sd in subdirs:
        if not os.path.isdir(sd):
            continue
            
        h_val, dt_val = extract_h_dt_from_folder(sd)
        if h_val is None or dt_val is None:
            continue
            
        pvtu_files = glob.glob(os.path.join(sd, "activation_time_*.pvtu"))
        if not pvtu_files:
            pvtu_files = glob.glob(os.path.join(sd, "activation_time_*.vtu"))
            
        if not pvtu_files:
            continue
            
        try:
            max_act = read_max_activation_time(pvtu_files)
            data_dict[(h_val, dt_val)] = max_act
            print(f"  [+] dt={dt_val:<6.4f} ms | h={h_val:<4.2f} mm -> Max Activation Time = {max_act:<6.4f} ms")
        except Exception as e:
            print(f"Error processing {sd}: {e}")

    if not data_dict:
        print("Error: No valid simulation data points extracted.")
        sys.exit(1)

    unique_h = np.array(sorted(list(set(k[0] for k in data_dict.keys()))))
    unique_dt = np.array(sorted(list(set(k[1] for k in data_dict.keys()))))
    
    # Grid indexing: X = Delta t, Y = Mesh Size h
    DT_grid, H_grid = np.meshgrid(unique_dt, unique_h, indexing='ij')
    Z_grid = np.zeros_like(DT_grid)

    for i, dt in enumerate(unique_dt):
        for j, h in enumerate(unique_h):
            Z_grid[i, j] = data_dict.get((h, dt), np.nan)

    dt_fine, h_fine, DT_fine, H_fine, Z_fine = interpolate_2d_grid(unique_dt, unique_h, Z_grid, num_x=100, num_y=100)

    # ----------------------------------------------------
    # Generate 3D Surface Plot (X = Delta t, Y = h)
    # ----------------------------------------------------
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection='3d')

    norm = mcolors.Normalize(vmin=Z_fine.min(), vmax=Z_fine.max())
    face_colors = cm.viridis(norm(Z_fine))

    surf = ax.plot_surface(DT_fine, H_fine, Z_fine, facecolors=face_colors, cmap='viridis', linewidth=0, antialiased=True, shade=False, alpha=0.85)
    sc = ax.scatter(DT_grid.ravel(), H_grid.ravel(), Z_grid.ravel(), c=Z_grid.ravel(), cmap='viridis', norm=norm, s=65, zorder=5, edgecolors='black', linewidth=0.8)

    for i in range(len(unique_dt)):
        for j in range(len(unique_h)):
            x, y, z = DT_grid[i, j], H_grid[i, j], Z_grid[i, j]
            ax.text(x, y, z + 0.6, f"{z:.2f}", fontsize=8.5, fontweight='bold', ha='center', zorder=6)

    ax.set_xlabel(r'Time Step $\Delta t$ [ms]', fontsize=11, labelpad=10)
    ax.set_ylabel('Mesh Size h [mm]', fontsize=11, labelpad=10)
    ax.set_zlabel('Activation Time [ms]', fontsize=11, labelpad=10)
    ax.invert_xaxis()
    ax.set_title(r'3D Convergence Sweep ($\Delta t$ on X-axis, h on Y-axis)' + '\nActivation Time vs. Time Step & Mesh Size', fontsize=13, pad=15, fontweight='bold')

    mappable = cm.ScalarMappable(norm=norm, cmap='viridis')
    mappable.set_array(Z_fine)
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label('Activation Time [ms]', fontsize=10)

    output_img = os.path.join(sweep_dir, "convergence_activation_time_3d.png")
    plt.tight_layout()
    plt.savefig(output_img, dpi=300, transparent=True)
    print(f"\n[Success] Swapped 3D graph saved successfully to: {output_img}")

if __name__ == "__main__":
    main()
