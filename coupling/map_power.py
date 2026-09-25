import os
import openmc
import numpy as np
import pandas as pd

# Force the script to run in its own directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def main():
  statepoint_path = '../neutronics/statepoint.100.h5'
  output_csv = 'power_map.csv'
  
  print(f"Loading {statepoint_path}...")
  sp = openmc.StatePoint(statepoint_path)
  
  # 1. Extract the power tally and its associated mesh
  tally = sp.get_tally(name='power_density')
  mesh_filter = tally.find_filter(openmc.MeshFilter)
  mesh = mesh_filter.mesh
  
  # Extract the raw tally mean values (eV per source particle)
  # OpenMC flattens mesh tallies such that x changes fastest, then y, then z.
  power_data = tally.mean.flatten()
  
  # 2. Reconstruct spatial coordinates (Centroids)
  nx, ny, nz = mesh.dimension
  ll = mesh.lower_left
  ur = mesh.upper_right
  
  # Calculate cell widths
  dx = (ur[0] - ll[0]) / nx
  dy = (ur[1] - ll[1]) / ny
  dz = (ur[2] - ll[2]) / nz
  
  # Generate 1D arrays for the cell centers
  x_centers = np.linspace(ll[0] + dx/2, ur[0] - dx/2, nx)
  y_centers = np.linspace(ll[1] + dy/2, ur[1] - dy/2, ny)
  z_centers = np.linspace(ll[2] + dz/2, ur[2] - dz/2, nz)
  
  # Create a 3D grid of coordinates matching OpenMC's flattening order (ij indexing)
  X, Y, Z = np.meshgrid(x_centers, y_centers, z_centers, indexing='ij')
  
  # 3. Normalize the power profile
  # For robust MOOSE coupling, it is best to provide a normalized shape factor
  # (where the sum of all points = 1.0) and multiply by the total reactor power 
  # natively inside the MOOSE input file.
  total_tally_sum = np.sum(power_data)
  if total_tally_sum > 0:
    normalized_power = power_data / total_tally_sum
  else:
    normalized_power = power_data
    print("Warning: Total power tally is zero.")
      
  # 4. Export to CSV for MOOSE (converting OpenMC cm to MOOSE meters)
  df = pd.DataFrame({
    'x': X.flatten() / 100.0,
    'y': Y.flatten() / 100.0,
    'z': Z.flatten() / 100.0,
    'power_fraction': normalized_power
  })
  
  df.to_csv(output_csv, index=False)
  print(f"Successfully mapped {len(df)} spatial points to {output_csv}")

if __name__ == "__main__":
  main()