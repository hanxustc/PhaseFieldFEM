#pragma once

#include "AuxKernel.h"
#include "qpNeighbors.h"

class RegGradAux : public AuxKernel
{
public:
  static InputParameters validParams();
  RegGradAux(const InputParameters & parameters);

protected:
  Real computeValue() override;
  
  const ADVariableGradient & _grad_theta;
  const qpNeighbors & _qp_neighbors;
  const unsigned int _component;
};