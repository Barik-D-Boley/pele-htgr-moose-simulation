# ==============================================================================
# CASE B -- Higher fidelity: power shape from the OpenMC neutronics solve
#
# First create openmc_power.txt (from the repository root):
#   python coupling/map_power.py
# Then run (from the moose/ folder):
#   mpiexec -n 4 combined-opt -i pele_core_sim.i case_openmc.i
# ==============================================================================

[Functions]
  [power_density]
    # Reads a 3D x-y-z grid of W/m^3 values and interpolates between points
    type = PiecewiseMultilinear
    data_file = openmc_power.txt
  []
[]

[Outputs]
  file_base = pele_openmc_out
[]
