import os
import subprocess as sp
from tqdm import tqdm

def generate_meshes():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    mesh_scripts_dir = os.path.join(script_dir, 'mesh_scripts')

    mesh_output_dir = os.path.join(script_dir, '..', 'mesh')
    os.makedirs(mesh_output_dir, exist_ok=True)

    geo_files = [f for f in os.listdir(mesh_scripts_dir) if f.endswith('.geo')]
    if not geo_files:
        print("No .geo files found in mesh_scripts.")
        return

    
    bar_format = "{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}"

    pbar = tqdm(
        geo_files,
        desc="Generating meshes",
        unit="file",
        bar_format=bar_format
    )

    for filename in pbar:
        if '2D' in filename:
            dimension = 2
        elif '3D' in filename:
            dimension = 3
        else:
            tqdm.write(f"Cannot determine dimension for {filename}, skipping.")
            continue
        
        pbar.set_description_str(f"Processing: {filename} ({dimension}D)")
        
        geo_path = os.path.join(mesh_scripts_dir, filename)
        output_name = filename.replace('.geo', '.msh')
        output_path = os.path.join(mesh_output_dir, output_name)

        cmd = f'gmsh -{dimension} -format msh "{geo_path}" -o "{output_path}"'
        result = sp.run(cmd, shell=True, stdout=sp.DEVNULL, stderr=sp.STDOUT)

        if result.returncode != 0:
            tqdm.write(f"Error generating mesh for {filename}.")
        

    tqdm.write("All meshes have been generated successfully.")

if __name__ == '__main__':
    generate_meshes()
