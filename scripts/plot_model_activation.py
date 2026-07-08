#!/usr/bin/env pvpython
"""
Script to generate clip plane visualizations of activation time for model sweeps.
Applies a ParaView Plane Clip centered at (10, 3.5, 1.5) with normal vector defined by
the cross product of (20, 3.5, 1.5) and (-20, 3.5, 0), colored using the 'Jet' colormap.
Camera is positioned starting from corner (0, 0, 3) looking towards the center (10, 3.5, 1.5).
All model outputs share a unified color scale spanning the global min and max activation times.
"""

import sys
import os
os.environ['DISPLAY'] = ''  # Force headless offscreen rendering to prevent window flashing
import glob
import re
import argparse
import numpy as np

try:
    from paraview.simple import *
    import paraview.servermanager as sm
    HAS_PARAVIEW = True
except ImportError:
    HAS_PARAVIEW = False

def find_latest_models_dir(base_dir="build"):
    candidates = glob.glob(os.path.join(base_dir, "*_models"))
    if not candidates:
        candidates = glob.glob("*_models")
    if not candidates:
        return None
    return max(candidates, key=os.path.getmtime)

def get_file_range(pvtu_path):
    """Extract scalar range for Activation Time from a VTU/PVTU file."""
    if pvtu_path.endswith('.pvtu'):
        reader = XMLPartitionedUnstructuredGridReader(FileName=[pvtu_path])
    else:
        reader = XMLUnstructuredGridReader(FileName=[pvtu_path])
        
    reader.UpdatePipeline()
    fetched = sm.Fetch(reader)
    arr = fetched.GetPointData().GetArray('Activation Time')
    rng = arr.GetRange() if arr else None
    Delete(reader)
    return rng

def render_model_clip(pvtu_path, output_png, min_val=None, max_val=None, show_legend=False, view_size=[1000, 1000]):
    """Render activation time on a clipped plane with Jet colormap and fixed color scale."""
    if not HAS_PARAVIEW:
        raise RuntimeError("ParaView python modules are not available.")

    # Load reader
    if pvtu_path.endswith('.pvtu'):
        reader = XMLPartitionedUnstructuredGridReader(FileName=[pvtu_path])
    else:
        reader = XMLUnstructuredGridReader(FileName=[pvtu_path])

    # Define Clip plane parameters
    # Origin: (10, 3.5, 1.5)
    origin = [10.0, 3.5, 1.5]
    
    # Normal: Cross product of (20, 3.5, 1.5) and (-20, 3.5, 0)
    v1 = np.array([20.0, 7, 3])
    v2 = np.array([-20.0, 7, 0.0])
    normal = list(np.cross(v1, v2))

    clip = Clip(Input=reader)
    clip.ClipType = 'Plane'
    clip.ClipType.Origin = origin
    clip.ClipType.Normal = normal

    # Create render view
    renderView = CreateView('RenderView')
    renderView.ViewSize = view_size
    renderView.Background = [0.15, 0.15, 0.15]  # Dark background for contrast
    renderView.OrientationAxesVisibility = 0

    display = Show(clip, renderView)
    display.Representation = 'Surface'
    ColorBy(display, ('POINTS', 'Activation Time'))

    # Apply Jet color preset and unify color scale across models
    lut = GetColorTransferFunction('ActivationTime')
    lut.ApplyPreset('Jet', True)
    if min_val is not None and max_val is not None:
        lut.RescaleTransferFunction(min_val, max_val)
    display.LookupTable = lut

    # Show colorbar scalar bar
    scalarBar = GetScalarBar(lut, renderView)
    scalarBar.Title = 'Activation Time [ms]'
    scalarBar.ComponentTitle = ''
    scalarBar.TitleColor = [0.0, 0.0, 0.0]
    scalarBar.LabelColor = [0.0, 0.0, 0.0]
    scalarBar.Visibility = 1 if show_legend else 0

    # Position camera at corner (0, 0, 3) looking towards domain center (10, 3.5, 1.5)
    renderView.CameraPosition = [-17.0525, -16.2422, 25.7423]
    renderView.CameraFocalPoint = origin
    renderView.CameraViewUp = [0.444443, 0.384503, 0.809091]

    Render()
    SaveScreenshot(output_png, renderView, ImageResolution=view_size, TransparentBackground=1)
    
    # Cleanup ParaView objects
    Delete(renderView)
    Delete(clip)
    Delete(reader)

def extract_activation_time_over_line(pvtu_path):
    """Extract activation time along the line (0,0,0) to (20,7,3)."""
    if not HAS_PARAVIEW:
        raise RuntimeError("ParaView python modules are not available.")

    # Load reader
    if pvtu_path.endswith('.pvtu'):
        reader = XMLPartitionedUnstructuredGridReader(FileName=[pvtu_path])
    else:
        reader = XMLUnstructuredGridReader(FileName=[pvtu_path])

    # Plot over line filter
    plot = PlotOverLine(Input=reader)
    plot.Point1 = [0.0, 0.0, 0.0]
    plot.Point2 = [20.0, 7.0, 3.0]

    plot.UpdatePipeline()

    # Fetch dataset to client
    fetched = sm.Fetch(plot)

    # Extract points and activation times using VTK-to-NumPy conversion
    from vtkmodules.util.numpy_support import vtk_to_numpy
    
    vtk_pts = fetched.GetPoints()
    if not vtk_pts:
        Delete(plot)
        Delete(reader)
        return None, None

    pts = vtk_to_numpy(vtk_pts.GetData())

    vtk_arr = fetched.GetPointData().GetArray('Activation Time')
    if not vtk_arr:
        Delete(plot)
        Delete(reader)
        return None, None

    act_time = vtk_to_numpy(vtk_arr)

    # Calculate 1D distance from point (0,0,0)
    dist = np.linalg.norm(pts, axis=1)

    # Sort arrays by distance to ensure clean line plotting
    sort_idx = np.argsort(dist)
    dist = dist[sort_idx]
    act_time = act_time[sort_idx]

    # Cleanup
    Delete(plot)
    Delete(reader)

    return dist, act_time

def plot_combined_activation_line(all_line_data, output_png):
    """Plot activation times along the line for all models on a shared figure."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Modern styling
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    fig, ax = plt.subplots(figsize=(10, 6))

    # Unified distinct colors for cardiac models
    model_colors = {
        'EPI': '#2b5c8f',    # Steel Blue
        'ENDO': '#e74c3c',   # Red/Coral
        'M': '#2ecc71',      # Green
        'PB': '#9b59b6',     # Yellow/Gold
        'TNNP': '#FFCD00'    # Purple
    }
    color_cycle = ['#1abc9c', '#34495e', '#d35400', '#7f8c8d']
    color_idx = 0

    for model_name, (dist, act_time) in sorted(all_line_data.items()):
        color = model_colors.get(model_name.upper())
        if not color:
            color = color_cycle[color_idx % len(color_cycle)]
            color_idx += 1

        # Plot only activated (non-NaN) values
        valid = ~np.isnan(act_time)
        if np.any(valid):
            ax.plot(dist[valid], act_time[valid], label=model_name, color=color, linewidth=2.5)

    ax.set_xlabel('Distance along line (0,0,0) -> (20,7,3) [mm]', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_ylabel('Activation Time [ms]', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_title('Activation Time Profile along Line (0,0,0) to (20,7,3)', fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(fontsize=10, frameon=True, facecolor='white', edgecolor='#bdc3c7')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#bdc3c7')
    ax.spines['bottom'].set_color('#bdc3c7')

    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.close(fig)
    print(f"  [+] Saved combined line plot to: {output_png}")

    # Plot 2: Standalone Zoomed Plot (x < 1.0)
    fig_zoom, ax_zoom = plt.subplots(figsize=(10, 6))
    for model_name, (dist, act_time) in sorted(all_line_data.items()):
        color = model_colors.get(model_name.upper())
        if not color:
            color = color_cycle[color_idx % len(color_cycle)]
            color_idx += 1
        valid = ~np.isnan(act_time)
        if np.any(valid):
            ax_zoom.plot(dist[valid], act_time[valid], label=model_name, color=color, linewidth=2.5)
            
    ax_zoom.set_xlabel('Distance along line (0,0,0) -> (20,7,3) [mm]', fontsize=12, fontweight='bold', labelpad=10)
    ax_zoom.set_ylabel('Activation Time [ms]', fontsize=12, fontweight='bold', labelpad=10)
    ax_zoom.set_title('Activation Time Profile (Zoomed x < 1.0 mm)\nLine (0,0,0) to (20,7,3)', fontsize=13, fontweight='bold', pad=15)
    ax_zoom.grid(True, linestyle='--', alpha=0.5)
    ax_zoom.legend(fontsize=10, frameon=True, facecolor='white', edgecolor='#bdc3c7')
    
    ax_zoom.set_xlim(0.0, 1.0)
    
    # Calculate limits for y-axis in the zoomed range
    y_vals_in_zoom = []
    for dist, act_time in all_line_data.values():
        valid_idx = (dist < 1.0) & (~np.isnan(act_time))
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
    output_png_zoom = output_png.replace(".png", "_zoom.png")
    plt.savefig(output_png_zoom, dpi=300)
    plt.close(fig_zoom)
    print(f"  [+] Saved zoomed line plot to: {output_png_zoom}")

    # Plot 3: Combined Main Plot with Zoomed Inset (x < 1.0)
    fig_inset, ax_main = plt.subplots(figsize=(10, 6))
    for model_name, (dist, act_time) in sorted(all_line_data.items()):
        color = model_colors.get(model_name.upper())
        if not color:
            color = color_cycle[color_idx % len(color_cycle)]
            color_idx += 1
        valid = ~np.isnan(act_time)
        if np.any(valid):
            ax_main.plot(dist[valid], act_time[valid], label=model_name, color=color, linewidth=2.5)
            
    ax_main.set_xlabel('Distance along line (0,0,0) -> (20,7,3) [mm]', fontsize=12, fontweight='bold', labelpad=10)
    ax_main.set_ylabel('Activation Time [ms]', fontsize=12, fontweight='bold', labelpad=10)
    ax_main.set_title('Activation Time Profile along Line (0,0,0) to (20,7,3)', fontsize=13, fontweight='bold', pad=15)
    ax_main.grid(True, linestyle='--', alpha=0.5)
    ax_main.legend(fontsize=10, frameon=True, facecolor='white', edgecolor='#bdc3c7')
    
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)
    ax_main.spines['left'].set_color('#bdc3c7')
    ax_main.spines['bottom'].set_color('#bdc3c7')
    
    # Position inset axes in the bottom-right corner [0.55, 0.15, 0.35, 0.35]
    ax_ins = ax_main.inset_axes([0.55, 0.15, 0.35, 0.35])
    for model_name, (dist, act_time) in sorted(all_line_data.items()):
        color = model_colors.get(model_name.upper())
        if not color:
            color = color_cycle[color_idx % len(color_cycle)]
            color_idx += 1
        valid = ~np.isnan(act_time)
        if np.any(valid):
            ax_ins.plot(dist[valid], act_time[valid], color=color, linewidth=2.0)
            
    ax_ins.set_xlim(0.0, 1.0)
    if y_vals_in_zoom:
        ymin, ymax = min(y_vals_in_zoom), max(y_vals_in_zoom)
        padding = (ymax - ymin) * 0.1 if ymax > ymin else 1.0
        ax_ins.set_ylim(ymin - padding, ymax + padding)
        
    ax_ins.grid(True, linestyle=':', alpha=0.6)
    
    try:
        ax_main.indicate_inset_zoom(ax_ins, edgecolor="black", alpha=0.3)
    except Exception:
        pass
        
    plt.tight_layout()
    output_png_combined = output_png.replace(".png", "_combined.png")
    plt.savefig(output_png_combined, dpi=300)
    plt.close(fig_inset)
    print(f"  [+] Saved combined line plot with inset to: {output_png_combined}")

def render_only_legend(pvtu_path, output_png, min_val=None, max_val=None):
    """Render only the colorbar legend directly inside ParaView by hiding the model."""
    if not HAS_PARAVIEW:
        raise RuntimeError("ParaView python modules are not available.")

    # Load reader
    if pvtu_path.endswith('.pvtu'):
        reader = XMLPartitionedUnstructuredGridReader(FileName=[pvtu_path])
    else:
        reader = XMLUnstructuredGridReader(FileName=[pvtu_path])

    # We need a clip/data source to create a scalar bar
    clip = Clip(Input=reader)
    clip.ClipType = 'Plane'

    renderView = CreateView('RenderView')
    renderView.ViewSize = [150, 400]
    renderView.Background = [0.15, 0.15, 0.15]  # Matches model image background
    renderView.OrientationAxesVisibility = 0

    display = Show(clip, renderView)
    display.Representation = 'Surface'
    ColorBy(display, ('POINTS', 'Activation Time'))

    lut = GetColorTransferFunction('ActivationTime')
    lut.ApplyPreset('Jet', True)
    if min_val is not None and max_val is not None:
        lut.RescaleTransferFunction(min_val, max_val)

    # Set representation opacity to 0.0 so the model is invisible,
    # but keep it shown so the colorbar legend renders its colors.
    display.Opacity = 0.0

    # Show and configure the colorbar legend
    scalarBar = GetScalarBar(lut, renderView)
    scalarBar.Title = 'Activation Time [ms]'
    scalarBar.ComponentTitle = ''
    scalarBar.TitleColor = [0.0, 0.0, 0.0]
    scalarBar.LabelColor = [0.0, 0.0, 0.0]
    scalarBar.TitleFontSize = 22
    scalarBar.LabelFontSize = 20
    scalarBar.Visibility = 1
    
    # Position colorbar nicely inside the 150x400 view
    scalarBar.Position = [0.1, 0.05]
    scalarBar.ScalarBarLength = 0.9

    Render()
    SaveScreenshot(output_png, renderView, ImageResolution=[150, 400], TransparentBackground=1)

    Delete(renderView)
    Delete(clip)
    Delete(reader)

def main():
    parser = argparse.ArgumentParser(description="Generate clipped activation time PNG plots for model sweeps.")
    parser.add_argument("models_dir", nargs="?", default=None, help="Path to the models sweep directory")
    args = parser.parse_args()

    models_dir = args.models_dir
    if not models_dir:
        models_dir = find_latest_models_dir()

    if not models_dir or not os.path.exists(models_dir):
        print(f"Error: Could not locate models sweep directory: {models_dir}")
        sys.exit(1)

    print(f"Processing model sweep directory: {models_dir}")

    subdirs = sorted(glob.glob(os.path.join(models_dir, "*")))
    tasks = []

    # 1. Collect master files only (.pvtu files preferred, skipping sub-partition pieces)
    for sd in subdirs:
        if not os.path.isdir(sd):
            continue
        pvtu_files = glob.glob(os.path.join(sd, "*.pvtu"))
        if not pvtu_files:
            all_vtu = glob.glob(os.path.join(sd, "*.vtu"))
            pvtu_files = [f for f in all_vtu if not re.search(r'\.[0-9]{2,}\.vtu$', f)]
            
        for pvtu in sorted(list(set(pvtu_files))):
            tasks.append((sd, pvtu))

    if not tasks:
        print("Warning: No VTU/PVTU master files found to process.")
        sys.exit(0)

    # 2. Compute Global Min and Max across all model datasets
    print("Computing global scalar range for unified colormap scale...")
    all_mins, all_maxs = [], []
    for sd, pvtu in tasks:
        rng = get_file_range(pvtu)
        if rng:
            all_mins.append(rng[0])
            all_maxs.append(rng[1])

    global_min = min(all_mins) if all_mins else None
    global_max = max(all_maxs) if all_maxs else None
    if global_min is not None and global_max is not None:
        print(f"  [+] Unified Color Scale Range: [{global_min:.4f}, {global_max:.4f}] ms")

    # 3. Render master screenshots with unified color scale
    # Extract the legend as a standalone image, then render all runs without legend as native squares
    legend_extracted = False
    for sd, pvtu in tasks:
        base_name = os.path.splitext(os.path.basename(pvtu))[0]
        output_png = os.path.join(models_dir, f"{base_name}_clip.png")
        
        # Extract colorbar legend directly using ParaView's hide feature if not done yet
        if not legend_extracted:
            colorbar_png = os.path.join(models_dir, "colorbar_legend.png")
            print(f"  [+] Natively rendering colorbar legend: {colorbar_png}")
            try:
                render_only_legend(pvtu, colorbar_png, min_val=global_min, max_val=global_max)
                legend_extracted = True
            except Exception as e:
                print(f"      -> Error rendering colorbar: {e}")
                
        print(f"  [+] Rendering (no legend, native square 1000x1000): {pvtu}")
        try:
            render_model_clip(pvtu, output_png, min_val=global_min, max_val=global_max, show_legend=False, view_size=[1000, 1000])
            print(f"      -> Saved PNG: {output_png}")
        except Exception as e:
            print(f"      -> Error rendering {pvtu}: {e}")

    # 4. Extract and plot activation time profiles along line (0,0,0) -> (20,7,3)
    print("Extracting activation time profiles along the line (0,0,0) -> (20,7,3)...")
    all_line_data = {}
    for sd, pvtu in tasks:
        model_name = os.path.basename(sd).split('_')[0]
        print(f"  [+] Extracting line profile for model: {model_name} ({pvtu})")
        try:
            dist, act_time = extract_activation_time_over_line(pvtu)
            if dist is not None and act_time is not None:
                all_line_data[model_name] = (dist, act_time)
        except Exception as e:
            print(f"      -> Error extracting line profile: {e}")

    if all_line_data:
        output_line_plot = os.path.join(models_dir, "activation_time_over_line.png")
        try:
            plot_combined_activation_line(all_line_data, output_line_plot)
        except Exception as e:
            print(f"      -> Error generating combined line plot: {e}")

if __name__ == "__main__":
    main()
