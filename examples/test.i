[GlobalParams]
  is_sharp = false
  const_M = true
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
  nx = 64
  ny = 64
  xmin = 0
  xmax = 10.0
  ymin = 0
  ymax = 10.0
  elem_type = QUAD4
  parallel_type = distributed
[]

[UserObjects]
  [./initial_theta]
    type = initialTheta
    theta_init = '-0.6 -0.45 -0.3 -0.15 0.0 0.15'
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
  [./grad_norm]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[Materials]
  # non-AD version: type = Functions
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
  [./gnorm]
    type = RegGradAux
    variable = grad_norm
    theta = theta
    qp_neighbors = qp_neighbors
    component = 2
    execute_on = 'initial timestep_end'
  [../]
[]

[Adaptivity]
  [./Markers]
    [./initial_uniform]
      type = UniformMarker
      mark = REFINE
    [../]
    [./gb_marker]
      type = ValueThresholdMarker
      variable = grad_norm
      refine = 0.1
      coarsen = 0.05
    [../]
  [../]
  marker = gb_marker
  initial_marker = initial_uniform
  max_h_level = 2
  initial_steps = 3
  interval = 500
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type -ksp_rtol -ksp_atol -ksp_max_it -snes_linesearch_type'
  petsc_options_value = 'lu       superlu_dist               1e-6      1e-9      100         bt'

  l_max_its = 100
  l_tol = 1e-6
  nl_max_its = 100
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-3

  dtmin = 1e-10
  start_time = 0.0
  end_time = 5.0
  dt = 1e-4
[]

[Postprocessors]
  [./total_energy]
    type = ADElementIntegralMaterialProperty
    mat_prop = energy
    execute_on = 'timestep_end FINAL'
  [../]
  [./dt]
    type = TimestepSize
  [../]
[]

[Outputs]
  [./console]
    type = Console
    max_rows = 10
    output_file = true
    time_step_interval = 100
  [../]
  [./exodus]
    type = Exodus
    time_step_interval = 1000
  [../]
  [./csv]
    type = CSV
    time_step_interval = 100
    execute_on = 'timestep_end FINAL'
  [../]
[]