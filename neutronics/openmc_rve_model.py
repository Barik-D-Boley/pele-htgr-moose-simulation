"""
OpenMC model of the Project Pele HTGR quarter-cell RVE.

Geometry (matches cad/htgr_rve_matrix.STEP and the Cubit mesh exactly):
  * Box: x = 0 -> 18.80 mm, y = 0 -> 9.40 mm, z = 0 -> 20.0 mm
  * Quarter COOLANT hole: radius 7.94 mm, centered on the corner (0, 0)
  * Quarter FUEL compact: radius 6.35 mm, centered on the corner (18.80, 9.40)
  * Everything else is H-451 graphite

OpenMC works in centimeters, so all dimensions below are in cm.
All six outer faces are reflective, which means "this cell repeats forever
in every direction" -- the standard way to model one cell of a big lattice.

Output: statepoint.<batches>.h5 containing k-eff and a 'fuel_power' tally
(fission energy deposited in each small piece of the fuel compact).
"""
import os

import numpy as np
import openmc

# Run inside the 'neutronics' folder so all outputs land here
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# ==============================================================================
# 0. DIMENSIONS (cm) -- keep in sync with the CAD / Cubit mesh
# ==============================================================================
CELL_X = 1.880      # RVE width  (x)
CELL_Y = 0.940      # RVE depth  (y)
HEIGHT = 2.000      # RVE height (z)

COOLANT_R = 0.794   # coolant hole radius, centered at (0, 0)
FUEL_R = 0.635      # fuel compact radius, centered at (CELL_X, CELL_Y)
FUEL_CENTER = (CELL_X, CELL_Y)

# ==============================================================================
# 1. MATERIALS
# ==============================================================================
# The fuel compact is thousands of tiny TRISO particles in a graphite binder.
# Here we "homogenize" it: smear the particles evenly through the compact.
# Each TRISO particle is a uranium kernel wrapped in 4 coating layers.
# (Dimensions follow the standard AGR-program TRISO design.)

# --- TRISO layer materials ---
kernel = openmc.Material(name='UO2 kernel (19.75% HALEU)')
kernel.add_element('U', 1.0, enrichment=19.75)
kernel.add_element('O', 2.0)
kernel.set_density('g/cm3', 10.4)

buffer = openmc.Material(name='Porous carbon buffer')
buffer.add_element('C', 1.0)
buffer.set_density('g/cm3', 1.05)

pyc = openmc.Material(name='Pyrolytic carbon (IPyC/OPyC)')
pyc.add_element('C', 1.0)
pyc.set_density('g/cm3', 1.90)

sic = openmc.Material(name='Silicon carbide')
sic.add_element('Si', 1.0)
sic.add_element('C', 1.0)
sic.set_density('g/cm3', 3.18)

compact_matrix = openmc.Material(name='Compact graphite binder')
compact_matrix.add_element('C', 1.0)
compact_matrix.set_density('g/cm3', 1.60)

# --- Outer radius of each TRISO layer (cm) ---
r_kernel = 0.02125            # 425 um diameter kernel
r_buffer = r_kernel + 0.0100  # 100 um buffer
r_ipyc = r_buffer + 0.0040    # 40 um inner PyC
r_sic = r_ipyc + 0.0035       # 35 um SiC
r_opyc = r_sic + 0.0040       # 40 um outer PyC

PACKING_FRACTION = 0.35       # fraction of compact volume filled by particles

# Volume fraction of each layer within the whole compact
radii = [0.0, r_kernel, r_buffer, r_ipyc, r_sic, r_opyc]
shell_fracs = [PACKING_FRACTION * (radii[i + 1]**3 - radii[i]**3) / r_opyc**3
               for i in range(5)]
layer_mats = [kernel, buffer, pyc, sic, pyc]

# PyC appears twice (inner + outer), so merge those fractions
mix = {}
for mat, frac in zip(layer_mats, shell_fracs):
    mix[mat] = mix.get(mat, 0.0) + frac
mix[compact_matrix] = 1.0 - PACKING_FRACTION

fuel = openmc.Material.mix_materials(list(mix.keys()), list(mix.values()), 'vo',
                                     name='Homogenized TRISO compact')
fuel.add_s_alpha_beta('c_Graphite')

graphite = openmc.Material(name='H-451 Graphite')
graphite.add_element('C', 1.0)
graphite.set_density('g/cm3', 1.74)
graphite.add_s_alpha_beta('c_Graphite')

helium = openmc.Material(name='Helium Coolant')
helium.add_element('He', 1.0)
helium.set_density('g/cm3', 0.0039)   # ~7 MPa, ~600 C

materials = openmc.Materials([fuel, graphite, helium])

# ==============================================================================
# 2. GEOMETRY
# ==============================================================================
x_min = openmc.XPlane(0.0, boundary_type='reflective')
x_max = openmc.XPlane(CELL_X, boundary_type='reflective')
y_min = openmc.YPlane(0.0, boundary_type='reflective')
y_max = openmc.YPlane(CELL_Y, boundary_type='reflective')
z_min = openmc.ZPlane(0.0, boundary_type='reflective')
z_max = openmc.ZPlane(HEIGHT, boundary_type='reflective')
box = +x_min & -x_max & +y_min & -y_max & +z_min & -z_max

coolant_cyl = openmc.ZCylinder(x0=0.0, y0=0.0, r=COOLANT_R)
fuel_cyl = openmc.ZCylinder(x0=FUEL_CENTER[0], y0=FUEL_CENTER[1], r=FUEL_R)

fuel_cell = openmc.Cell(name='fuel', fill=fuel, region=-fuel_cyl & box)
coolant_cell = openmc.Cell(name='coolant', fill=helium, region=-coolant_cyl & box)
graphite_cell = openmc.Cell(name='graphite', fill=graphite,
                            region=+fuel_cyl & +coolant_cyl & box)

geometry = openmc.Geometry([fuel_cell, coolant_cell, graphite_cell])

# ==============================================================================
# 3. SETTINGS
# ==============================================================================
settings = openmc.Settings()
settings.batches = 150      # total batches
settings.inactive = 50      # batches thrown away while the fission source settles
settings.particles = 20000  # neutrons per batch (more = less statistical noise)

# Start neutrons uniformly in the fuel compact region
source_box = openmc.stats.Box(
    [FUEL_CENTER[0] - FUEL_R, FUEL_CENTER[1] - FUEL_R, 0.0],
    [FUEL_CENTER[0], FUEL_CENTER[1], HEIGHT])
settings.source = openmc.IndependentSource(space=source_box,
                                           constraints={'fissionable': True})

# ==============================================================================
# 4. TALLIES (for MOOSE coupling)
# ==============================================================================
# A cylindrical mesh wrapped around the fuel axis, covering ONLY the fuel.
# This way no tally bin is half-fuel/half-graphite, so the power shape near
# the compact edge is captured accurately.
#   r:   10 rings from the axis out to the compact edge
#   phi: 6 slices across the quarter (180 deg -> 270 deg points into the cell)
#   z:   5 axial slices
fuel_mesh = openmc.CylindricalMesh(
    r_grid=np.linspace(0.0, FUEL_R, 11),
    phi_grid=np.linspace(np.pi, 1.5 * np.pi, 7),
    z_grid=np.linspace(0.0, HEIGHT, 6),
    origin=(FUEL_CENTER[0], FUEL_CENTER[1], 0.0),
    name='fuel_mesh')

power_tally = openmc.Tally(name='fuel_power')
power_tally.filters = [openmc.MeshFilter(fuel_mesh)]
power_tally.scores = ['kappa-fission']   # recoverable fission energy

tallies = openmc.Tallies([power_tally])

# ==============================================================================
# 5. RUN
# ==============================================================================
model = openmc.Model(geometry=geometry, materials=materials,
                     settings=settings, tallies=tallies)

if __name__ == '__main__':
    statepoint_path = model.run()
    with openmc.StatePoint(statepoint_path) as sp:
        print(f'\nk-eff = {sp.keff}')
    print(f'Results written to {statepoint_path}')
