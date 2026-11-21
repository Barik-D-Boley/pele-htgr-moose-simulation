# ==============================================================================
# Project Pele: HTGR RVE (quarter cell)
# Coupled Steady-State Thermal-Mechanical Analysis -- SHARED BASE INPUT
#
# This file holds everything the two cases have in common. It is NOT run on
# its own. Combine it with ONE case file, which supplies the heat source:
#
#   Simple assumption (uniform heat in the fuel):
#     mpiexec -n 4 combined-opt -i pele_core_sim.i case_uniform.i
#
#   Higher fidelity (OpenMC power shape):
#     mpiexec -n 4 combined-opt -i pele_core_sim.i case_openmc.i
#
# Both cases put the SAME total power into the fuel; only the shape differs.
# Each case file defines two functions:
#   power_density        -- fission heating in the fuel (W/m^3)
#   coolant_temperature  -- helium temperature along the channel, T(z) in K
# Helium flows DOWN the channel: 250 C at the top (z = 0.8 m), 700 C at the
# bottom (z = 0). How fast it heats up depends on the case's power shape.
# ==============================================================================

[Mesh]
  # SI-scaled (meters) Exodus mesh from meshes/create_mesh.jou
  # Blocks:  matrix_graphite, fuel_compact
  [fmg]
    type = FileMeshGenerator
    file = '../meshes/htgr_rve_mesh.e'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [temperature]
    order = FIRST
    family = LAGRANGE
    initial_condition = 750 # K, starting guess near the average coolant temperature
  []
[]

[Physics]
  # Automatically adds displacement variables and stress-divergence kernels
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = true
        eigenstrain_names = thermal_expansion
        generate_output = 'vonmises_stress hydrostatic_stress'
        temperature = temperature
      []
    []
  []
[]

[AuxVariables]
  # Height of each node, used to report WHERE the peak temperature occurs
  [z_position]
    order = FIRST
    family = LAGRANGE
  []
  # Copy of the heat source, so it can be plotted next to the results
  [power_q]
    order = CONSTANT
    family = MONOMIAL
    block = 'fuel_compact'
  []
[]

[AuxKernels]
  [z_position]
    type = ParsedAux
    variable = z_position
    expression = 'z'
    use_xyzt = true
  []
  [power_q]
    type = FunctionAux
    variable = power_q
    function = power_density
    block = 'fuel_compact'
  []
[]

[Kernels]
  # Steady heat equation:  -div(k grad T) = Q
  [heat_conduction]
    type = HeatConduction
    variable = temperature
  []

  # Volumetric fission heating, ONLY inside the fuel compact.
  # The function 'power_density' (W/m^3) is defined in the case file.
  [fission_heating]
    type = HeatSource
    variable = temperature
    block = 'fuel_compact'
    function = power_density
  []
[]

[BCs]
  # Helium cooling on the coolant channel wall. The helium temperature
  # varies with height (material property 'T_coolant', defined below).
  [convective_cooling]
    type = ConvectiveHeatFluxBC
    variable = temperature
    boundary = 'coolant_wall'
    T_infinity = T_coolant
    # Heat transfer coefficient CONSISTENT with the helium flow rate.
    # The flow is whatever carries the fuel power away with a 250 -> 700 C
    # rise: only ~1.3 g/s per channel, so the flow is barely turbulent
    # (Reynolds ~2800) and h is low. coupling/map_power.py prints this
    # estimate ("Sanity check: ... estimated h") -- if you change the power
    # or the inlet/outlet temperatures, update this number to match.
    # (Earlier placeholder was 2000, which would need ~10x more flow.)
    heat_transfer_coefficient = 150 # W/m^2-K
  []

  # All other outer faces: no heat flow (default), because neighboring cells
  # are identical (symmetry).

  # --- Symmetry roller constraints (stop the part from floating away) ---
  [symm_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'symm_x_min'
    value = 0
  []
  [symm_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'symm_y_min'
    value = 0
  []
  [bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'bottom_face'
    value = 0
  []
[]

[Materials]
  # Helium temperature on the coolant wall, from the case's coolant function
  [coolant]
    type = GenericFunctionMaterial
    block = 'matrix_graphite'
    prop_names = 'T_coolant'
    prop_values = 'coolant_temperature'
  []

  # ---------------------------------------------------------------------------
  # H-451 nuclear graphite (matrix block)
  # ---------------------------------------------------------------------------
  [graphite_conductivity]
    type = ParsedMaterial
    block = 'matrix_graphite'
    property_name = 'thermal_conductivity'
    coupled_variables = 'temperature'
    expression = '11000 / temperature' # W/m-K, simple fit for H-451
  []
  [graphite_elasticity]
    type = ComputeIsotropicElasticityTensor
    block = 'matrix_graphite'
    youngs_modulus = 9.8e9
    poissons_ratio = 0.13
  []
  [graphite_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = 'matrix_graphite'
    temperature = temperature
    stress_free_temperature = 300
    thermal_expansion_coeff = 4.0e-6
    eigenstrain_name = thermal_expansion
  []

  # ---------------------------------------------------------------------------
  # TRISO fuel compact (fuel block)
  # PLACEHOLDER values typical of TRISO compacts -- refine from literature.
  # ---------------------------------------------------------------------------
  [fuel_conductivity]
    type = GenericConstantMaterial
    block = 'fuel_compact'
    prop_names = 'thermal_conductivity'
    prop_values = '15.0' # W/m-K
  []
  [fuel_elasticity]
    type = ComputeIsotropicElasticityTensor
    block = 'fuel_compact'
    youngs_modulus = 15.0e9
    poissons_ratio = 0.15
  []
  [fuel_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = 'fuel_compact'
    temperature = temperature
    stress_free_temperature = 300
    thermal_expansion_coeff = 5.0e-6
    eigenstrain_name = thermal_expansion
  []

  # Stress calculation for both blocks
  [linear_stress]
    type = ComputeLinearElasticStress
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_rel_tol = 1e-8
[]

[Postprocessors]
  # --- Results to compare between the two cases ---
  [max_T_fuel]
    type = NodalExtremeValue
    variable = temperature
    block = 'fuel_compact'
  []
  [max_T_graphite]
    type = NodalExtremeValue
    variable = temperature
    block = 'matrix_graphite'
  []
  [max_vonmises_fuel]
    type = ElementExtremeValue
    variable = vonmises_stress
    block = 'fuel_compact'
  []
  [max_vonmises_graphite]
    type = ElementExtremeValue
    variable = vonmises_stress
    block = 'matrix_graphite'
  []

  # --- Where along the column is the fuel hottest? ---
  # Reports the height (z, in m) of the node where the fuel temperature is
  # highest: it looks up 'z_position' at the node where 'temperature' peaks.
  [z_of_max_T_fuel]
    type = NodalExtremeValue
    variable = z_position
    proxy_variable = temperature
    block = 'fuel_compact'
    value_type = max
  []

  # --- Energy conservation check: power_in should equal power_out ---
  [power_in]
    type = FunctionElementIntegral
    function = power_density
    block = 'fuel_compact'
  []
  [power_out]
    type = SideDiffusiveFluxIntegral
    variable = temperature
    diffusivity = thermal_conductivity
    boundary = 'coolant_wall'
  []

  # --- Compute time (seconds) for the cost comparison ---
  [run_time]
    type = PerfGraphData
    section_name = 'Root'
    data_type = TOTAL
  []
[]

[VectorPostprocessors]
  # Profiles along the full column height, read by
  # post-processing/compare_cases.py. Each line runs bottom (z = 0) to top.
  # 161 points = every 5 mm.

  # Near the fuel compact axis (hottest fuel)
  [fuel_centerline]
    type = LineValueSampler
    variable = 'temperature vonmises_stress power_q'
    start_point = '0.0186 0.0092 0.0'
    end_point = '0.0186 0.0092 0.8'
    num_points = 161
    sort_by = z
  []
  # In the graphite just outside the fuel wall, on the side facing the
  # coolant hole (fuel/graphite interface stresses)
  [graphite_near_fuel]
    type = LineValueSampler
    variable = 'temperature vonmises_stress'
    start_point = '0.0122 0.0092 0.0'
    end_point = '0.0122 0.0092 0.8'
    num_points = 161
    sort_by = z
  []
  # In the graphite just outside the coolant wall, on the 45-degree line
  [graphite_near_coolant]
    type = LineValueSampler
    variable = 'temperature vonmises_stress'
    start_point = '0.0058 0.0058 0.0'
    end_point = '0.0058 0.0058 0.8'
    num_points = 161
    sort_by = z
  []
[]

[Outputs]
  # file_base is set in each case file so the results don't overwrite
  exodus = true
  csv = true
  perf_graph = true
  print_linear_residuals = false
[]
