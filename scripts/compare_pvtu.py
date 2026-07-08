#!/usr/bin/env python3
"""
Compare PVTu files from different codes/runs.

This script compares solution fields from multiple PVTu files and computes
various similarity metrics including L2 norm, max absolute difference, 
and relative errors.

Usage:
    python compare_pvtu.py file1.pvtu file2.pvtu [--field solution] [--tolerance 1e-6]
    python compare_pvtu.py --dir1 output1/ --dir2 output2/ [--pattern "*.pvtu"]
"""

import argparse
import sys
from pathlib import Path
from typing import Tuple, Dict, List
import numpy as np


def load_pvtu_field(path: Path, field_name: str = "solution"):
    """Load a specific field from a PVTu or VTU file."""
    try:
        import pyvista as pv
    except ImportError:
        print("Error: pyvista is required. Install with: pip install pyvista")
        sys.exit(1)
    
    mesh = pv.read(str(path))
    
    # Try point data first, then cell data
    if field_name in mesh.point_data:
        return mesh.point_data[field_name], mesh, 'point'
    elif field_name in mesh.cell_data:
        return mesh.cell_data[field_name], mesh, 'cell'
    else:
        available_point = list(mesh.point_data.keys())
        available_cell = list(mesh.cell_data.keys())
        raise ValueError(
            f"Field '{field_name}' not found.\n"
            f"Available point data: {available_point}\n"
            f"Available cell data: {available_cell}"
        )


def compute_metrics(arr1: np.ndarray, arr2: np.ndarray) -> Dict[str, float]:
    """Compute similarity metrics between two arrays."""
    
    # Ensure arrays have compatible shapes
    if arr1.shape != arr2.shape:
        raise ValueError(f"Shape mismatch: {arr1.shape} vs {arr2.shape}")
    
    diff = arr1 - arr2
    abs_diff = np.abs(diff)
    
    metrics = {
        'max_abs_diff': np.max(abs_diff),
        'mean_abs_diff': np.mean(abs_diff),
        'l2_norm_diff': np.linalg.norm(diff),
        'l2_norm_ref': np.linalg.norm(arr1),
        'relative_l2_error': np.linalg.norm(diff) / (np.linalg.norm(arr1) + 1e-16),
        'max_abs_value_1': np.max(np.abs(arr1)),
        'max_abs_value_2': np.max(np.abs(arr2)),
        'mean_value_1': np.mean(arr1),
        'mean_value_2': np.mean(arr2),
        'relative_mean_diff': abs(np.mean(arr1) - np.mean(arr2)) / (abs(np.mean(arr1)) + 1e-16),
    }
    
    # Pointwise relative error (avoid division by zero)
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_error = abs_diff / (np.abs(arr1) + 1e-16)
        metrics['max_relative_error'] = np.max(rel_error)
        metrics['mean_relative_error'] = np.mean(rel_error)
    
    return metrics


def compare_two_files(file1: Path, file2: Path, field_name: str = "solution", 
                     tolerance: float = 1e-6, verbose: bool = True) -> bool:
    """Compare two PVTu files and return True if they are similar within tolerance."""
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"Comparing:")
        print(f"  File 1: {file1}")
        print(f"  File 2: {file2}")
        print(f"{'='*70}")
    
    try:
        data1, mesh1, loc1 = load_pvtu_field(file1, field_name)
        data2, mesh2, loc2 = load_pvtu_field(file2, field_name)
    except Exception as e:
        print(f"Error loading files: {e}")
        return False
    
    if verbose:
        print(f"\nField '{field_name}' loaded from {loc1} data")
        print(f"  Array shape: {data1.shape} vs {data2.shape}")
        print(f"  Data type: {data1.dtype} vs {data2.dtype}")
    
    try:
        metrics = compute_metrics(data1, data2)
    except Exception as e:
        print(f"Error computing metrics: {e}")
        return False
    
    # Print metrics
    if verbose:
        print(f"\n{'Metric':<30} {'Value':>15}")
        print(f"{'-'*45}")
        print(f"{'Max absolute difference':<30} {metrics['max_abs_diff']:>15.6e}")
        print(f"{'Mean absolute difference':<30} {metrics['mean_abs_diff']:>15.6e}")
        print(f"{'L2 norm of difference':<30} {metrics['l2_norm_diff']:>15.6e}")
        print(f"{'L2 norm of reference':<30} {metrics['l2_norm_ref']:>15.6e}")
        print(f"{'Relative L2 error':<30} {metrics['relative_l2_error']:>15.6e}")
        print(f"{'Max relative error':<30} {metrics['max_relative_error']:>15.6e}")
        print(f"{'Mean relative error':<30} {metrics['mean_relative_error']:>15.6e}")
        print(f"\n{'Mean value (file 1)':<30} {metrics['mean_value_1']:>15.6e}")
        print(f"{'Mean value (file 2)':<30} {metrics['mean_value_2']:>15.6e}")
        print(f"{'Relative mean difference':<30} {metrics['relative_mean_diff']:>15.6e}")
    
    # Determine if files are similar
    is_similar = metrics['relative_l2_error'] < tolerance
    
    if verbose:
        print(f"\n{'='*70}")
        if is_similar:
            print(f"✓ FILES ARE SIMILAR (relative L2 error < {tolerance})")
        else:
            print(f"✗ FILES ARE DIFFERENT (relative L2 error >= {tolerance})")
        print(f"{'='*70}\n")
    
    return is_similar


def compare_directories(dir1: Path, dir2: Path, pattern: str = "*.pvtu",
                       field_name: str = "solution", tolerance: float = 1e-6) -> Dict:
    """Compare all matching PVTu files in two directories."""
    
    files1 = sorted(dir1.glob(pattern))
    files2 = sorted(dir2.glob(pattern))
    
    if not files1:
        print(f"No files matching '{pattern}' found in {dir1}")
        return {}
    
    if not files2:
        print(f"No files matching '{pattern}' found in {dir2}")
        return {}
    
    print(f"\nFound {len(files1)} files in directory 1")
    print(f"Found {len(files2)} files in directory 2")
    
    # Match files by name
    files1_dict = {f.name: f for f in files1}
    files2_dict = {f.name: f for f in files2}
    
    common_names = set(files1_dict.keys()) & set(files2_dict.keys())
    only_in_1 = set(files1_dict.keys()) - set(files2_dict.keys())
    only_in_2 = set(files2_dict.keys()) - set(files1_dict.keys())
    
    if only_in_1:
        print(f"\nFiles only in directory 1: {sorted(only_in_1)}")
    if only_in_2:
        print(f"Files only in directory 2: {sorted(only_in_2)}")
    
    if not common_names:
        print("\nNo common files found between directories!")
        return {}
    
    print(f"\nComparing {len(common_names)} common files...\n")
    
    results = {}
    similar_count = 0
    
    for name in sorted(common_names):
        file1 = files1_dict[name]
        file2 = files2_dict[name]
        
        is_similar = compare_two_files(file1, file2, field_name, tolerance, verbose=False)
        results[name] = is_similar
        
        if is_similar:
            similar_count += 1
            status = "✓ SIMILAR"
        else:
            status = "✗ DIFFERENT"
        
        print(f"{status:15} {name}")
    
    print(f"\n{'='*70}")
    print(f"Summary: {similar_count}/{len(common_names)} files are similar (tolerance={tolerance})")
    print(f"{'='*70}\n")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Compare PVTu files from different codes/runs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare two specific files
  python compare_pvtu.py file1.pvtu file2.pvtu
  
  # Compare with custom field name
  python compare_pvtu.py file1.pvtu file2.pvtu --field temperature
  
  # Compare all files in two directories
  python compare_pvtu.py --dir1 output1/ --dir2 output2/
  
  # Compare with custom tolerance
  python compare_pvtu.py file1.pvtu file2.pvtu --tolerance 1e-3
  
  # Compare specific pattern
  python compare_pvtu.py --dir1 run1/ --dir2 run2/ --pattern "*_solution_*.pvtu"
        """
    )
    
    # Two modes: individual files or directories
    parser.add_argument('file1', nargs='?', type=Path,
                       help='First PVTu file to compare')
    parser.add_argument('file2', nargs='?', type=Path,
                       help='Second PVTu file to compare')
    parser.add_argument('--dir1', type=Path,
                       help='First directory containing PVTu files')
    parser.add_argument('--dir2', type=Path,
                       help='Second directory containing PVTu files')
    parser.add_argument('--pattern', default='*.pvtu',
                       help='Glob pattern for files (default: *.pvtu)')
    parser.add_argument('--field', default='solution',
                       help='Field name to compare (default: solution)')
    parser.add_argument('--tolerance', type=float, default=1e-6,
                       help='Relative L2 error tolerance (default: 1e-6)')
    parser.add_argument('--list-fields', action='store_true',
                       help='List available fields in the first file and exit')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.dir1 or args.dir2:
        if not (args.dir1 and args.dir2):
            parser.error("Both --dir1 and --dir2 must be provided for directory comparison")
        if args.file1 or args.file2:
            parser.error("Cannot specify both individual files and directories")
        mode = 'directory'
    elif args.file1 and args.file2:
        mode = 'files'
    else:
        parser.error("Must specify either two files or two directories to compare")
    
    # List fields mode
    if args.list_fields:
        if mode == 'files':
            test_file = args.file1
        else:
            test_file = next(args.dir1.glob(args.pattern), None)
            if not test_file:
                print(f"No files matching '{args.pattern}' found in {args.dir1}")
                sys.exit(1)
        
        try:
            import pyvista as pv
            mesh = pv.read(str(test_file))
            print(f"\nAvailable fields in {test_file}:")
            print(f"\nPoint data fields:")
            for key in mesh.point_data.keys():
                print(f"  - {key}")
            print(f"\nCell data fields:")
            for key in mesh.cell_data.keys():
                print(f"  - {key}")
        except Exception as e:
            print(f"Error reading file: {e}")
            sys.exit(1)
        sys.exit(0)
    
    # Compare files or directories
    if mode == 'files':
        if not args.file1.exists():
            print(f"Error: File not found: {args.file1}")
            sys.exit(1)
        if not args.file2.exists():
            print(f"Error: File not found: {args.file2}")
            sys.exit(1)
        
        is_similar = compare_two_files(args.file1, args.file2, 
                                      args.field, args.tolerance)
        sys.exit(0 if is_similar else 1)
    
    else:  # directory mode
        if not args.dir1.is_dir():
            print(f"Error: Directory not found: {args.dir1}")
            sys.exit(1)
        if not args.dir2.is_dir():
            print(f"Error: Directory not found: {args.dir2}")
            sys.exit(1)
        
        results = compare_directories(args.dir1, args.dir2, args.pattern,
                                     args.field, args.tolerance)
        
        # Exit code: 0 if all common files are similar, 1 otherwise
        all_similar = all(results.values()) if results else False
        sys.exit(0 if all_similar else 1)


if __name__ == '__main__':
    main()
