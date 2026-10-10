#include "ADGrad.h"
#include "libmesh/utility.h"
#include "libmesh/libmesh_common.h"

registerMooseObject("PhaseFieldApp", ADGrad);

InputParameters
ADGrad::validParams()
{
  InputParameters params = ADKernelGrad::validParams();
  params.addRequiredCoupledVar("theta", "The global orientation variable");
  params.addRequiredCoupledVar("thetaP", "Nonlocal theta plus");
  params.addRequiredCoupledVar("thetaM", "Nonlocal theta minus");

  params.addRequiredParam<bool>("is_sharp", "Sharp cut-off?");
  params.addRequiredParam<Real>("alpha", "Constant for additional smoothness");
  params.addRequiredParam<Real>("beta", "Constant for bulk constraint");
  params.addRequiredParam<Real>("u0", "Lower bound for GBs");
  
  return params;
}

ADGrad::ADGrad(const InputParameters & params)
  : ADKernelGrad(params),
    _theta(adCoupledValue("theta")),
    _grad_theta(adCoupledGradient("theta")),
    _thetaP(coupledValue("thetaP")),
    _thetaM(coupledValue("thetaM")),

    _is_sharp(getParam<bool>("is_sharp")),
    _alpha(getParam<Real>("alpha")),
    _beta(getParam<Real>("beta")),
    _u0(getParam<Real>("u0")),

    _B(getADMaterialProperty<Real>("B")),
    _R(getADMaterialProperty<Real>("R")),
    _w(getADMaterialProperty<Real>("w")),
    _c(getADMaterialProperty<Real>("c")),
    _dB(getADMaterialProperty<RealVectorValue>("dB")),
    _dw(getADMaterialProperty<RealVectorValue>("dw")),
    _dc(getADMaterialProperty<Real>("dc")),
    _psi(getADMaterialProperty<Real>("psi")),
    _M(getADMaterialProperty<Real>("M"))
{}

ADRealVectorValue
ADGrad::precomputeQpResidual()
{
  const Real eps = 1e-6;
  Real mask = 1.0;

  if (_is_sharp) {
    if (_grad_theta[_qp].norm() <= _u0) return 0.0;
  } else {
    Real grad = MetaPhysicL::raw_value(_grad_theta[_qp].norm());
    mask = 0.5 * (1.0 + std::tanh((grad - _u0)/eps));
  }
  
  Real dtheta = _thetaP[_qp] - _thetaM[_qp];
  Real coeff = 1.0 / (pow(dtheta, 2) + eps);

  ADRealVectorValue df(0.0, 0.0, 0.0);
  df = _M[_qp] * (_B[_qp]*coeff*(2*_w[_qp]*_grad_theta[_qp] + _dw[_qp]*(_grad_theta[_qp]*_grad_theta[_qp])) + _dB[_qp]*coeff*_w[_qp]*(_grad_theta[_qp]*_grad_theta[_qp])
      + coeff * _alpha * _R[_qp] * (2*(1 - _w[_qp])*_grad_theta[_qp] - _dw[_qp]*(_grad_theta[_qp]*_grad_theta[_qp]))
      - _beta * _R[_qp] * _dw[_qp] * _c[_qp]);

  return mask * df;
}