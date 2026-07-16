#include "RegGradAux.h"

registerMooseObject("PhaseFieldApp", RegGradAux);

InputParameters
RegGradAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("theta", "Orientation variable");
  params.addRequiredParam<UserObjectName>("qp_neighbors", "qpNeighbors userobject");
  params.addRequiredParam<unsigned int>("component", "0=x, 1=y");
  return params;
}

RegGradAux::RegGradAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _grad_theta(adCoupledGradient("theta")),
    _qp_neighbors(getUserObject<qpNeighbors>("qp_neighbors")),
    _component(getParam<unsigned int>("component"))
{
}

Real
RegGradAux::computeValue()
{
  const auto g = _qp_neighbors.useSmooth()
    ? _qp_neighbors.getSmoothGrad(_current_elem->id(), _qp)
    : MetaPhysicL::raw_value(_grad_theta[_qp]);
  if (_component == 0)
    return g(0);
  else if (_component == 1)
    return g(1);

  return 0.0;
}