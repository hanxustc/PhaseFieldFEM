#pragma once
#include "ADKernel.h"
#include "ADKernelGrad.h"

class ADGrad : public ADKernelGrad
{
public:
  static InputParameters validParams();
  ADGrad(const InputParameters & params);

protected:
  const ADVariableValue & _theta;
  const ADVariableGradient & _grad_theta;
  const VariableValue & _thetaP;
  const VariableValue & _thetaM;

  const bool _is_sharp;
  const Real _alpha;
  const Real _beta;
  const Real _u0;
  
  const ADMaterialProperty<Real> & _B;
  const ADMaterialProperty<Real> & _R;
  const ADMaterialProperty<Real> & _w;
  const ADMaterialProperty<Real> & _c;
  const ADMaterialProperty<RealVectorValue> & _dB;
  const ADMaterialProperty<RealVectorValue> & _dw;
  const ADMaterialProperty<Real> & _dc;
  const ADMaterialProperty<Real> & _psi;
  const ADMaterialProperty<Real> & _M;
  
  virtual ADRealVectorValue precomputeQpResidual() override;
};