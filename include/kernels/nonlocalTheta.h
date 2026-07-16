#pragma once

#include "AuxKernel.h"
#include "initialTheta.h"
#include "qpNeighbors.h"

class nonlocalTheta : public AuxKernel
{
public:
  static InputParameters validParams();
  nonlocalTheta(const InputParameters & params);

protected:
  virtual Real computeValue() override;

  const ADVariableValue & _theta;
  const ADVariableGradient & _grad_theta;

  const qpNeighbors & _qp_neighbors;
  const initialTheta & _initial_theta;
  const std::vector<Real> & _theta_init;
  
  Real thetaP;
  Real thetaM;
  const Real _u0;
  const int _every;
  const bool _is_plus;
};