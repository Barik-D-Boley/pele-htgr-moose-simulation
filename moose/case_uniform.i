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

  # Uniform power -> the helium heats up in a straight line.
  # Flow is DOWN: inlet at the top (z = H), outlet at the bottom (z = 0).
  # T_in, T_out must match --t-in, --t-out in coupling/map_power.py.
  [coolant_temperature]
    type = ParsedFunction
    expression = 'T_in + (T_out - T_in) * (H - z) / H'
    symbol_names = 'T_in    T_out   H'
    symbol_values = '523.15 973.15 0.8' # K, K, m
  []
[]

[Outputs]
  file_base = pele_uniform_out
[]
