[GlobalParams]
  is_sharp = false
  const_M = false
  M0 = 1.0
  alpha = 0.01
  beta = 1.0
  u0 = 0.0001
  u1_0 = 32.0
  every = 50
  dmax = 1.0
  dmin = 0.2
  use_smooth = false
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 80
  ny = 80
  nz = 0
  xmin = 0
  xmax = 10.0
  ymin = 0
  ymax = 10.0
  zmin = 0
  zmax = 0
  elem_type = QUAD4
  uniform_refine = 0
  
  parallel_type = distributed
[]

[UserObjects]
  [./initial_theta]
    type = initialTheta
    theta_init = '-0.05 0.0 0.05 0.15 0.3 0.6'
    execute_on = 'initial'
  [../]

  [./qp_neighbors]
    type = qpNeighbors
    theta = theta
    initial_theta = initial_theta
    execute_on = 'initial timestep_begin'
  [../]
[]

[Variables]
  [./theta]
    family = LAGRANGE
    order = FIRST
  [../]
[]

[AuxVariables]
  [./thetaP]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./thetaM]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./grad_x]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./grad_y]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[Materials]
  [./functions]
    type = ADFunctions
    theta = theta
    thetaP = thetaP
    thetaM = thetaM
    B0 = 1.0
    n = 2
    m = 4
    eps_m = 0.0
    K = 0.0
  [../]
[]

[Functions]
  [./theta_from_file]
    type = PiecewiseMultilinear
    data_file = poly.data
  [../]
[]

[ICs]
  [./ic]
    type = FunctionIC
    variable = theta
    function = theta_from_file
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      variable = theta
      auto_direction = 'x y'
    [../]
  [../]
[]

[Kernels]
  [./time]
    type = ADTimeDerivative
    variable = theta
  [../]
  
  [./gradient]
    type = ADGrad
    variable = theta
    theta = theta
    thetaP = thetaP
    thetaM = thetaM
  [../]

  [./bulk]
    type = ADBulk
    variable = theta
    theta = theta
  [../]
[]

[AuxKernels]
  [./thetaP]
    type = nonlocalTheta
    variable = thetaP
    initial_theta = initial_theta
    theta = theta
    qp_neighbors = qp_neighbors
    is_plus = true
    execute_on = 'initial timestep_end'
  [../]
  [./thetaM]
    type = nonlocalTheta
    variable = thetaM
    initial_theta = initial_theta
    theta = theta
    qp_neighbors = qp_neighbors
    is_plus = false
    execute_on = 'initial timestep_end'
  [../]
[]

#[Adaptivity]
#  [./Indicators]
#    [./grad_theta]
#      type = GradientJumpIndicator
#      variable = theta
#    [../]
#  [../]

#  [./Markers]
#    [./gb_marker]
#      type = ErrorFractionMarker
#      indicator = grad_theta
#      refine = 0.2
#      coarsen = 0.05
#    [../]
#  [../]

#  marker = gb_marker
#  max_h_level = 3
#  initial_steps = 3
#  interval = 500
#[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'NEWTON'

  petsc_options_iname = '-pc_type -ksp_rtol -ksp_atol -ksp_max_it -snes_linesearch_type'
  petsc_options_value = 'lu       1e-6      1e-9      100         bt'

  l_max_its = 100 # Max number of linear iterations
  l_tol = 1e-6 # Relative tolerance for linear solves
  nl_max_its = 100 # Max number of nonlinear iterations
  nl_abs_tol = 1e-6 # Absolute tolerance for nonlienar solves
  nl_rel_tol = 1e-3 # Relative tolerance for nonlienar solves

  dtmin = 1e-6
  dtmax = 1e-4
  start_time = 0.0
  end_time = 20.0

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-4
    optimal_iterations = 9
    iteration_window = 3
    growth_factor = 1.2
    cutback_factor = 0.5
  [../]
[]

[Postprocessors]
  [./total_energy]
    type = ADElementIntegralMaterialProperty
    mat_prop = energy
    execute_on = 'timestep_end FINAL'
  [../]
  [./energy_change]
    type = ChangeOverTimePostprocessor
    postprocessor = total_energy
    execute_on = 'timestep_end FINAL'
  [../]
[]

[UserObjects]
  [./terminator]
    type = Terminator
    expression = 'abs(energy_change) < 1e-9'
    message = 'Energy converged'
    execute_on = 'timestep_end'
  [../]
[]

[Outputs]
  [./checkpoint]
    type = Checkpoint
    num_files = 20
    time_step_interval = 10000
    file_base = /oscar/data/avandewa/xhan39/Projects/pf_fem/test_out
  [../]

  [./console]
    type = Console
    max_rows = 10
    output_file = true
    time_step_interval = 100
  [../]
  
  [./exodus]
    type = Exodus
    time_step_interval = 100
  [../]
  
  [./csv]
    type = CSV
    time_step_interval = 100
    execute_on = 'timestep_end FINAL'
  [../]
[]