import os
import openmc

# Force the script to run in its own directory (the 'neutronics' folder)
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# ==============================================================================
# 1. MATERIALS
# ==============================================================================
haleu = openmc.Material(name='HALEU TRISO')
haleu.add_nuclide('U235', 0.1975)
haleu.add_nuclide('U238', 0.8025)
haleu.add_element('O', 2.0)
haleu.set_density('g/cm3', 10.5)

graphite = openmc.Material(name='H-451 Graphite')
graphite.add_element('C', 1.0)
graphite.set_density('g/cm3', 1.74)
graphite.add_s_alpha_beta('c_Graphite')

helium = openmc.Material(name='Helium Coolant')
helium.add_element('He', 1.0)
helium.set_density('g/cm3', 0.0015)

materials = openmc.Materials([haleu, graphite, helium])
materials.export_to_xml()

# ==============================================================================
# 2. GEOMETRY (18.80 mm x 9.40 mm RVE)
# ==============================================================================
pitch_x = 1.88  # cm
pitch_y = 0.94  # cm
height = 2.0    # cm

min_x = openmc.XPlane(x0=-pitch_x/2, boundary_type='reflective')
max_x = openmc.XPlane(x0=pitch_x/2, boundary_type='reflective')
min_y = openmc.YPlane(y0=-pitch_y/2, boundary_type='reflective')
max_y = openmc.YPlane(y0=pitch_y/2, boundary_type='reflective')
min_z = openmc.ZPlane(z0=-height/2, boundary_type='reflective')
max_z = openmc.ZPlane(z0=height/2, boundary_type='reflective')

fuel_cyl = openmc.ZCylinder(r=0.3)

fuel_cell = openmc.Cell(fill=haleu, region=-fuel_cyl & +min_z & -max_z)
mod_cell = openmc.Cell(fill=graphite, region=+fuel_cyl & +min_x & -max_x & +min_y & -max_y & +min_z & -max_z)

universe = openmc.Universe(cells=[fuel_cell, mod_cell])
geometry = openmc.Geometry(universe)
geometry.export_to_xml()

# ==============================================================================
# 3. SETTINGS
# ==============================================================================
settings = openmc.Settings()
settings.batches = 100
settings.inactive = 20
settings.particles = 1000

# Initial source distribution
space = openmc.stats.Box([-pitch_x/2, -pitch_y/2, -height/2], 
                         [pitch_x/2, pitch_y/2, height/2])
settings.source = openmc.IndependentSource(space=space)
settings.export_to_xml()

# ==============================================================================
# 4. TALLIES (For MOOSE Coupling)
# ==============================================================================
# Create a 3D mesh for spatial power mapping
mesh = openmc.RegularMesh()
mesh.dimension = [10, 10, 10]
mesh.lower_left = [-pitch_x/2, -pitch_y/2, -height/2]
mesh.upper_right = [pitch_x/2, pitch_y/2, height/2]

power_tally = openmc.Tally(name='power_density')
power_tally.filters = [openmc.MeshFilter(mesh)]
power_tally.scores = ['kappa-fission']

tallies = openmc.Tallies([power_tally])
tallies.export_to_xml()

# ==============================================================================
# 5. RUN
# ==============================================================================
openmc.run()