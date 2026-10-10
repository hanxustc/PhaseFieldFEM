#pragma once
#include "ADKernel.h"
#include "ADFunctions.h"

class ADBulk : public ADKernel
{
public:
  static InputParameters validParams();
  ADBulk(const InputParameters & params);

protected:
  const bool _is_sharp;
  const Real _beta;
  const Real _u0;
  
  const ADVariableValue & _theta;
  const ADVariableGradient & _grad_theta;
  const ADMaterialProperty<Real> & _B;
  const ADMaterialProperty<Real> & _R;
  const ADMaterialProperty<Real> & _w;
  const ADMaterialProperty<Real> & _c;
  const ADMaterialProperty<Real> & _dc;
  const ADMaterialProperty<Real> & _M;
  
  virtual ADReal computeQpResidual() override;
};