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
    initial_condition = 600 # K, starting guess near coolant temperature
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
  # Helium cooling on the coolant channel wall
  [convective_cooling]
    type = ConvectiveHeatFluxBC
    variable = temperature
    boundary = 'coolant_wall'
    T_infinity = 600              # K, bulk helium temperature
    heat_transfer_coefficient = 2000 # W/m^2-K
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

[Outputs]
  # file_base is set in each case file so the results don't overwrite
  exodus = true
  csv = true
  perf_graph = true
  print_linear_residuals = false
[]
