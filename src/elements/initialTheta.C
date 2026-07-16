#include "initialTheta.h"

registerMooseObject("PhaseFieldApp", initialTheta);

InputParameters
initialTheta::validParams()
{
  InputParameters params = NodalUserObject::validParams();
  params.addParam<std::vector<Real>>("theta_init", "List of initial theta values");
  return params;
}

initialTheta::initialTheta(const InputParameters & params)
  : NodalUserObject(params),
  _theta_init(getParam<std::vector<Real>>("theta_init"))
{}

void
initialTheta::initialize()
{}

void
initialTheta::execute()
{}

void
initialTheta::threadJoin(const UserObject & other)
{}

void
initialTheta::finalize()
{
  if (processor_id() == 0) {
    std::cout << "initialTheta: found " << _theta_init.size() << " distinct theta values:" << std::endl;
    for (const auto & t : _theta_init)
      std::cout << "  " << t;
    std::cout << std::endl;
  }
}

const std::vector<Real> &
initialTheta::getThetaInit() const
{
  return _theta_init;
}