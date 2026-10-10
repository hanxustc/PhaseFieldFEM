#include "Bulk.h"

registerMooseObject("PhaseFieldApp", Bulk);

InputParameters
Bulk::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("theta", "The global orientation variable");

  params.addRequiredParam<bool>("is_sharp", "Sharp cut-off?");
  params.addRequiredParam<Real>("beta", "Constant for bulk constraint");
  params.addRequiredParam<Real>("u0", "Lower bound for GBs");
  return params;
}

Bulk::Bulk(const InputParameters & params)
  : Kernel(params),
    _is_sharp(getParam<bool>("is_sharp")),
    _beta(getParam<Real>("beta")),
    _u0(getParam<Real>("u0")),

    _theta(coupledValue("theta")),
    _grad_theta(coupledGradient("theta")),

    _B(getMaterialProperty<Real>("B")),
    _R(getMaterialProperty<Real>("R")),
    _w(getMaterialProperty<Real>("w")),
    _c(getMaterialProperty<Real>("c")),
    _dc(getMaterialProperty<Real>("dc")),
    _M(getMaterialProperty<Real>("M")),

    _dw(getMaterialProperty<RealVectorValue>("dw")),
    _d2c(getMaterialProperty<Real>("d2c"))
{}

Real
Bulk::computeQpResidual()
{
  Real eps = 1e-6;
  Real mask = 1.0;

  if (_is_sharp) {
    if (_grad_theta[_qp].norm() <= _u0) return 0.0;
  } else {
    mask = 0.5 * (1.0 + std::tanh((_grad_theta[_qp].norm() - _u0)/eps));
  }

  return _test[_i][_qp] * _M[_qp] * _beta * _R[_qp] * (1 - _w[_qp]) * mask * _dc[_qp];
}

Real
Bulk::computeQpJacobian()
{
  Real eps = 1e-6;
  Real mask = 1.0;

  if (_is_sharp) {
    if (_grad_theta[_qp].norm() <= _u0) return 0.0;
  } else {
    mask = 0.5 * (1.0 + std::tanh((_grad_theta[_qp].norm() - _u0)/eps));
  }
  
  return _test[_i][_qp] * _M[_qp] * _beta * mask * _R[_qp] * ((1 - _w[_qp]) * _d2c[_qp] * _phi[_j][_qp] - (_dw[_qp] * _grad_phi[_j][_qp]) * _dc[_qp]);
}
