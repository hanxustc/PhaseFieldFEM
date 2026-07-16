#include "ADFunctions.h"

registerMooseObject("PhaseFieldApp", ADFunctions);

InputParameters
ADFunctions::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addRequiredCoupledVar("theta", "The global orientation variable");
  params.addRequiredCoupledVar("thetaP", "Nonlocal theta plus");
  params.addRequiredCoupledVar("thetaM", "Nonlocal theta minus");

  params.addRequiredParam<bool>("const_M", "Is mobility const?");
  params.addRequiredParam<Real>("M0", "Mobility");
  params.addRequiredParam<Real>("alpha", "Constant for additional smoothness");
  params.addRequiredParam<Real>("beta", "Constant for bulk constraint");
  params.addRequiredParam<Real>("u0", "Lower bound for GBs");
  params.addRequiredParam<Real>("u1_0", "Upper bound for GBs");
  params.addRequiredParam<Real>("B0", "Constant for anisotropic GB function");
  params.addRequiredParam<int>("n", "Symmetry for misorientation");
  params.addRequiredParam<int>("m", "Symmetry for Inclination");
  params.addRequiredParam<Real>("eps_m", "Amplitude of anisotropy");
  params.addRequiredParam<Real>("K", "Initial inclination");

  return params;
}

ADFunctions::ADFunctions(const InputParameters & parameters)
  : ADMaterial(parameters),
    _theta(adCoupledValue("theta")),
    _grad_theta(adCoupledGradient("theta")),
    _thetaP(coupledValue("thetaP")),
    _thetaM(coupledValue("thetaM")),

    _B(declareADProperty<Real>("B")),
    _R(declareADProperty<Real>("R")),
    _w(declareADProperty<Real>("w")),
    _c(declareADProperty<Real>("c")),
    _dB(declareADProperty<RealVectorValue>("dB")),
    _dw(declareADProperty<RealVectorValue>("dw")),
    _dc(declareADProperty<Real>("dc")),
    _psi(declareADProperty<Real>("psi")),
    _M(declareADProperty<Real>("M")),
    _energy(declareADProperty<Real>("energy")),

    _const_M(getParam<bool>("const_M")),
    _M0(getParam<Real>("M0")),
    _alpha(getParam<Real>("alpha")),
    _beta(getParam<Real>("beta")),
    _u0(getParam<Real>("u0")),
    _u1_0(getParam<Real>("u1_0")),
    _B0(getParam<Real>("B0")),
    _n(getParam<int>("n")),
    _m(getParam<int>("m")),
    _eps_m(getParam<Real>("eps_m")),
    _K(getParam<Real>("K"))
{
  _U(0,0) = 0;
  _U(0,1) = 1;
  _U(1,0) = -1;
  _U(1,1) = 0;
}

void
ADFunctions::computeQpProperties()
{
  const Real eps = 1e-6;
  
  Real dtheta = _thetaP[_qp] - _thetaM[_qp];
  ADReal u = 1.0;
  if (dtheta > eps)
    u = std::max(0.0, std::min(1.0, (_theta[_qp] - _thetaM[_qp]) / dtheta));

  // ψ = arctan(grad_theta_y / grad_theta_x)
  _psi[_qp] = 0.0;
  if (_grad_theta[_qp].norm() > _u0)
    _psi[_qp] = std::atan2(_grad_theta[_qp](1), _grad_theta[_qp](0));

  // Mobility
  if (_const_M) _M[_qp] = _M0;
  else {
    const Real c0 = 0.2;
    const Real half_period = M_PI / _n;
    const Real quarter_period = M_PI / (2.0 * _n);
    Real dt = std::fmod(dtheta, half_period);
    if (dt < quarter_period) {
      _M[_qp] = _M0 * (dt / c0) * (dt / c0) * (1.0 - std::exp(-dt / c0));
    } else {
      _M[_qp] = _M0 * ((half_period - dt) / c0) * ((half_period - dt) / c0) * (1.0 - std::exp(-(half_period - dt) / c0));
    }
  }

  // GB functions
  const Real theta_star = 0.5 * (_thetaP[_qp] + _thetaM[_qp]) + _K;
  _R[_qp] = std::abs(std::sin(_n * dtheta));
  _B[_qp] = _B0 * std::abs(std::sin(_n * dtheta)) * (1.0 + _eps_m*std::sin(_m * (theta_star - _psi[_qp])));

  _dB[_qp] = _B0 * std::abs(std::sin(_n * dtheta)) * _m * _eps_m * std::cos(_m * (theta_star - _psi[_qp])) 
                 / (_grad_theta[_qp]*_grad_theta[_qp] + eps) * _U * _grad_theta[_qp];
  
  // Weighting function
  const Real _u1 = _u1_0 * (dtheta > eps ? dtheta : eps);
  if (_grad_theta[_qp].norm() < _u0)
    _w[_qp] = 0.0;
  else if (_grad_theta[_qp].norm() > _u1)
    _w[_qp] = 1.0;
  else
    _w[_qp] = pow(_grad_theta[_qp].norm() - _u0, 2) * pow(_grad_theta[_qp].norm() - 2*_u1 + _u0, 2) / pow(_u1 - _u0, 4);

  _dw[_qp] = 0.0;
  if (_grad_theta[_qp].norm() >= _u0 && _grad_theta[_qp].norm() <= _u1)
    _dw[_qp] = 4 * (_grad_theta[_qp].norm() - _u0) * (_grad_theta[_qp].norm() - 2*_u1 + _u0) * (_grad_theta[_qp].norm() - _u1)
                 / pow(_u1 - _u0, 4) * _grad_theta[_qp]
                 / _grad_theta[_qp].norm();

  // Double-well potential
  _c[_qp] = pow(u, 2) * pow(1 - u, 2);
  _dc[_qp] = 2 * u * (1 - 3*u + 2*pow(u, 2)) / (dtheta + eps);

  // Energy density
  Real coeff = 0.0;
  if (dtheta > eps)
    coeff = 1 / dtheta;

  _energy[_qp] = _B[_qp]*_w[_qp]*_grad_theta[_qp]*_grad_theta[_qp]*pow(coeff, 2) + _R[_qp]*(1 - _w[_qp])*(_alpha*_grad_theta[_qp]*_grad_theta[_qp]*pow(coeff, 2) + _beta*_c[_qp]);
}