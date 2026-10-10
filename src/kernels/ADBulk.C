#include "ADBulk.h"

registerMooseObject("PhaseFieldApp", ADBulk);

InputParameters
ADBulk::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addRequiredCoupledVar("theta", "The global orientation variable");

  params.addRequiredParam<bool>("is_sharp", "Sharp cut-off?");
  params.addRequiredParam<Real>("beta", "Constant for bulk constraint");
  params.addRequiredParam<Real>("u0", "Lower bound for GBs");
  return params;
}

ADBulk::ADBulk(const InputParameters & params)
  : ADKernel(params),
    _theta(adCoupledValue("theta")),
    _grad_theta(adCoupledGradient("theta")),

    _is_sharp(getParam<bool>("is_sharp")),
    _beta(getParam<Real>("beta")),
    _u0(getParam<Real>("u0")),

    _B(getADMaterialProperty<Real>("B")),
    _R(getADMaterialProperty<Real>("R")),
    _w(getADMaterialProperty<Real>("w")),
    _c(getADMaterialProperty<Real>("c")),
    _dc(getADMaterialProperty<Real>("dc")),
    _M(getADMaterialProperty<Real>("M"))
{}

ADReal
ADBulk::computeQpResidual()
{
  Real eps = 1e-6;
  Real mask = 1.0;

  if (_is_sharp) {
    if (_grad_theta[_qp].norm() <= _u0) return 0.0;
  } else {
    Real grad = MetaPhysicL::raw_value(_grad_theta[_qp].norm());
    mask = 0.5 * (1.0 + std::tanh((grad - _u0)/eps));
  }

  return _test[_i][_qp] * _M[_qp] * _beta * _R[_qp] * (1 - _w[_qp]) * mask * _dc[_qp];
}