# ==============================================================================
# CASE A -- Simplified assumption: uniform heat generation in the fuel
#
# Run (from the moose/ folder):
#   mpiexec -n 4 combined-opt -i pele_core_sim.i case_uniform.i
# ==============================================================================

[Functions]
  [power_density]
    type = ConstantFunction
    value = 3.0e7 # W/m^3 -- must match --q-avg used in coupling/map_power.py
  []
[]

[Outputs]
  file_base = pele_uniform_out
[]
