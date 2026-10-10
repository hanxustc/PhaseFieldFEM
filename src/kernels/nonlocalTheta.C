#include "nonlocalTheta.h"

registerMooseObject("PhaseFieldApp", nonlocalTheta);

InputParameters
nonlocalTheta::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<UserObjectName>("initial_theta", "initialTheta nodal object");
  params.addRequiredParam<UserObjectName>("qp_neighbors", "qpNeighbors user object");
  params.addRequiredCoupledVar("theta", "The global orientation variable");
  
  params.addRequiredParam<Real>("u0", "Lower bound for GBs");
  params.addRequiredParam<int>("every", "Timestep interval for updating P/M neighbors");
  params.addRequiredParam<bool>("is_plus", "true -> compute thetaP (max), false -> thetaM (min)");

  return params;
}

nonlocalTheta::nonlocalTheta(const InputParameters & parameters)
  : AuxKernel(parameters),
    _theta(adCoupledValue("theta")),
    _grad_theta(adCoupledGradient("theta")),

    _initial_theta(getUserObject<initialTheta>("initial_theta")),
    _qp_neighbors(getUserObject<qpNeighbors>("qp_neighbors")),
    _theta_init(_initial_theta.getThetaInit()),

    _u0(getParam<Real>("u0")),
    _every(getParam<int>("every")),
    _is_plus(getParam<bool>("is_plus"))
{}

template <typename T>
static inline std::pair<T, T>
weight(const libMesh::VectorValue<T> & norm_grad, const libMesh::VectorValue<T> & u)
{
  const T t  = norm_grad * u;
  const T t0 = 0.8;
  const T t1 = 1.0;
  const T n  = 0.5;

  T Wp = 0.0;
  T Wm = 0.0;

  if (t > t0 && t < t1) {
    const T a = std::exp(-1.0 / std::pow(t - t0, n));
    const T b = std::exp(-1.0 / std::pow(t1 - t, n));
    Wp = a / (a + b);
  } else if (t <= t0)
    Wp = 0.0;
  else
    Wp = 1.0;

  if (t > -t1 && t < -t0) {
    const T a = std::exp(-1.0 / std::pow(t + t1, n));
    const T b = std::exp(-1.0 / std::pow(-t0 - t, n));
    Wm = 1.0 - a/(a + b);
  } else if (t <= -t1)
    Wm = 1.0;
  else
    Wm = 0.0;

  return {Wp, Wm};
}

Real
closestTheta(Real theta, const std::vector<Real> & theta_init)
{
  Real best = theta_init.front();
  Real min_dist = std::abs(theta - best);

  for (const auto & t : theta_init) {
    const Real dd = std::abs(theta - t);
    if (dd < min_dist) {
      min_dist = dd;
      best = t;
    }
  }
  return best;
}

Real
nonlocalTheta::computeValue()
{
  if (_theta_init.size() > 1) {
    const libMesh::Point & p0 = _q_point[_qp];
    const Real theta0 = MetaPhysicL::raw_value(_theta[_qp]);
    const Real grad_norm = _qp_neighbors.useSmooth() ? _qp_neighbors.getSmoothGradNorm(_current_elem->id(), _qp) : MetaPhysicL::raw_value(_grad_theta[_qp].norm());

    // if (!_qp_neighbors.neighborsRebuilt())
    if (_t_step % _every != 0)
      return _u[_qp];

    if (grad_norm <= 100*_u0)
      return theta0;
    
    const auto & all_qps = _qp_neighbors.getAllQps();
    const auto & neighP = _qp_neighbors.getPlusNeighbors(_current_elem->id(), _qp);
    const auto & neighM = _qp_neighbors.getMinusNeighbors(_current_elem->id(), _qp);
    const std::size_t n_layers = std::max(neighP.size(), neighM.size());
    bool plus_in_bulk  = false;
    bool minus_in_bulk = false;
    std::size_t last_layer = 0;

    for (std::size_t layer = 0; layer < n_layers; layer++) {
      last_layer  = layer;

      if (layer < neighP.size()) {
        const auto & nbP = all_qps[neighP[layer]];
        const Real nbP_norm = _qp_neighbors.useSmooth() ? _qp_neighbors.getSmoothGradNorm(nbP.elem_id, nbP.qp_id) : nbP.grad_theta.norm();
        if (nbP_norm <= _u0)
          plus_in_bulk = true;
      }

      if (layer < neighM.size()) {
        const auto & nbM = all_qps[neighM[layer]];
        const Real nbM_norm = _qp_neighbors.useSmooth() ? _qp_neighbors.getSmoothGradNorm(nbM.elem_id, nbM.qp_id) : nbM.grad_theta.norm();
        if (nbM_norm <= _u0)
          minus_in_bulk = true;
      }

      if (plus_in_bulk && minus_in_bulk)
        break;
    }

    if (neighP.empty() && neighM.empty()) {
      thetaP = theta0;
      thetaM = theta0;
    } else {
      thetaP = std::max(all_qps[neighP[last_layer]].theta, all_qps[neighM[last_layer]].theta);
      thetaM = std::min(all_qps[neighP[last_layer]].theta, all_qps[neighM[last_layer]].theta);
    }

    if (theta0 < thetaM)
      thetaM = theta0;
    if (theta0 > thetaP)
      thetaP = theta0;

    thetaP = closestTheta(thetaP, _theta_init);
    thetaM = closestTheta(thetaM, _theta_init);
  } else if (_theta_init.size() <= 1) {
    thetaP = *std::max_element(_theta_init.begin(), _theta_init.end());
    thetaM = *std::min_element(_theta_init.begin(), _theta_init.end());
  }
  
  /*
  libMesh::VectorValue<Real> norm_gradtheta(0, 0, 0);
  norm_gradtheta = MetaPhysicL::raw_value(_grad_theta[_qp]) / MetaPhysicL::raw_value(_grad_theta[_qp].norm());

  Real Wp_sum = 0.0;
  Real Wm_sum = 0.0;
  for (const auto & nb : neigh_qps) {
    auto v = _qp_neighbors.periodDirect(p0, nb.point);
    const Real normv = v.norm();
    if (normv > 1e-16)
      v /= normv;

    const auto [Wp_temp, Wm_temp] = weight(norm_gradtheta, v);
    Wp_sum += Wp_temp;
    Wm_sum += Wm_temp;
  }

  const Real inv_Wp = (Wp_sum > 1e-16) ? 1.0 / Wp_sum : 0.0;
  const Real inv_Wm = (Wm_sum > 1e-16) ? 1.0 / Wm_sum : 0.0;

  Real theta1 = 0.0;
  Real theta2 = 0.0;
  
  for (const auto & nb : neigh_qps) {
    auto v = _qp_neighbors.periodDirect(p0, nb.point);
    const Real normv = v.norm();
    if (normv > 1e-16)
      v /= normv;
    
    const auto [Wp_raw, Wm_raw] = weight(norm_gradtheta, v);
    const Real Wp = Wp_raw * inv_Wp;
    const Real Wm = Wm_raw * inv_Wm;
    
    theta1 += Wp * MetaPhysicL::raw_value(nb.theta);
    theta2 += Wm * MetaPhysicL::raw_value(nb.theta);
  }

  thetaP = std::max(theta1, theta2);
  thetaM = std::min(theta1, theta2);
  */

  return _is_plus ? thetaP : thetaM;
}