#include "Grad.h"
#include "libmesh/utility.h"
#include "libmesh/libmesh_common.h"

registerMooseObject("PhaseFieldApp", Grad);

InputParameters
Grad::validParams()
{
  InputParameters params = KernelGrad::validParams();
  params.addRequiredCoupledVar("theta", "The global orientation variable");
  params.addRequiredCoupledVar("thetaP", "Nonlocal theta plus");
  params.addRequiredCoupledVar("thetaM", "Nonlocal theta minus");

  params.addRequiredParam<bool>("is_sharp", "Sharp cut-off?");
  params.addRequiredParam<Real>("alpha", "Constant for additional smoothness");
  params.addRequiredParam<Real>("beta", "Constant for bulk constraint");
  params.addRequiredParam<Real>("u0", "Lower bound for GBs");

  return params;
}

Grad::Grad(const InputParameters & params)
  : KernelGrad(params),
    _theta(coupledValue("theta")),
    _grad_theta(coupledGradient("theta")),
    _thetaP(coupledValue("thetaP")),
    _thetaM(coupledValue("thetaM")),

    _is_sharp(getParam<bool>("is_sharp")),
    _alpha(getParam<Real>("alpha")),
    _beta(getParam<Real>("beta")),
    _u0(getParam<Real>("u0")),

    _B(getMaterialProperty<Real>("B")),
    _R(getMaterialProperty<Real>("R")),
    _w(getMaterialProperty<Real>("w")),
    _c(getMaterialProperty<Real>("c")),
    _dB(getMaterialProperty<RealVectorValue>("dB")),
    _dw(getMaterialProperty<RealVectorValue>("dw")),
    _dc(getMaterialProperty<Real>("dc")),
    _psi(getMaterialProperty<Real>("psi")),
    _M(getMaterialProperty<Real>("M")),

    _d2B(getMaterialProperty<RankTwoTensor>("d2B")),
    _d2w(getMaterialProperty<RankTwoTensor>("d2w"))
{
  _I(0,0) = 1;
  _I(0,1) = 0;
  _I(1,0) = 0;
  _I(1,1) = 1;
}

RealGradient
Grad::precomputeQpResidual()
{
  const Real eps = 1e-6;
  Real mask = 1.0;

  if (_is_sharp) {
    if (_grad_theta[_qp].norm() <= _u0) return 0.0;
  } else {
    mask = 0.5 * (1.0 + std::tanh((_grad_theta[_qp].norm() - _u0)/eps));
  }

  Real dtheta = _thetaP[_qp] - _thetaM[_qp];
  Real coeff = 1.0 / (pow(dtheta, 2) + eps);

  RealGradient df(0.0, 0.0, 0.0);
  df = _M[_qp] * (_B[_qp]*coeff*(2*_w[_qp]*_grad_theta[_qp] + _dw[_qp]*(_grad_theta[_qp]*_grad_theta[_qp])) + _dB[_qp]*coeff*_w[_qp]*(_grad_theta[_qp]*_grad_theta[_qp])
      + coeff * _alpha * _R[_qp] * (2*(1 - _w[_qp])*_grad_theta[_qp] - _dw[_qp]*(_grad_theta[_qp]*_grad_theta[_qp]))
      - _beta * _R[_qp] * _dw[_qp] * _c[_qp]);

  return mask * df;
}

RealGradient
Grad::precomputeQpJacobian()
{
  const Real eps = 1e-6;
  Real mask = 1.0;

  if (_is_sharp) {
    if (_grad_theta[_qp].norm() <= _u0) return 0.0;
  } else {
    mask = 0.5 * (1.0 + std::tanh((_grad_theta[_qp].norm() - _u0)/eps));
  }

  Real dtheta = _thetaP[_qp] - _thetaM[_qp];
  Real coeff = 1.0 / (pow(dtheta, 2) + eps);
  RankTwoTensor d2f;
  d2f = _M[_qp] * (coeff * (RankTwoTensor::outerProduct(2*_w[_qp]*_grad_theta[_qp] + _dw[_qp]*(_grad_theta[_qp]*_grad_theta[_qp]), _dB[_qp])
                           + _B[_qp]*(2*RankTwoTensor::outerProduct(_grad_theta[_qp], _dw[_qp]) + 2*_w[_qp]*_I + _d2w[_qp]*(_grad_theta[_qp]*_grad_theta[_qp]) + 2*RankTwoTensor::outerProduct(_dw[_qp], _grad_theta[_qp]))
                           + _d2B[_qp]*_w[_qp]*(_grad_theta[_qp]*_grad_theta[_qp])
                           + RankTwoTensor::outerProduct(_dB[_qp], _dw[_qp])*(_grad_theta[_qp]*_grad_theta[_qp])
                           + 2*_w[_qp]*RankTwoTensor::outerProduct(_dB[_qp], _grad_theta[_qp]))
                  + coeff * _alpha * _R[_qp] * (2*(1 - _w[_qp])*_I - 2*RankTwoTensor::outerProduct(_grad_theta[_qp], _dw[_qp]) - _d2w[_qp]*(_grad_theta[_qp]*_grad_theta[_qp]) - 2*RankTwoTensor::outerProduct(_dw[_qp], _grad_theta[_qp]))
                  - _beta * _R[_qp] * _d2w[_qp] * _c[_qp]);

  return mask * (d2f * _grad_phi[_j][_qp] - _M[_qp] * _beta * _R[_qp] * _dc[_qp] * _dw[_qp] * _phi[_j][_qp]);
}
