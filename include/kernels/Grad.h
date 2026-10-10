#pragma once
#include "KernelGrad.h"
#include "RankTwoTensor.h"

class Grad : public KernelGrad
{
public:
  static InputParameters validParams();
  Grad(const InputParameters & params);

protected:
  const VariableValue & _theta;
  const VariableGradient & _grad_theta;
  const VariableValue & _thetaP;
  const VariableValue & _thetaM;

  const bool _is_sharp;
  const Real _alpha;
  const Real _beta;
  const Real _u0;

  const MaterialProperty<Real> & _B;
  const MaterialProperty<Real> & _R;
  const MaterialProperty<Real> & _w;
  const MaterialProperty<Real> & _c;
  const MaterialProperty<RealVectorValue> & _dB;
  const MaterialProperty<RealVectorValue> & _dw;
  const MaterialProperty<Real> & _dc;
  const MaterialProperty<Real> & _psi;
  const MaterialProperty<Real> & _M;
  const MaterialProperty<RankTwoTensor> & _d2B;
  const MaterialProperty<RankTwoTensor> & _d2w;

  virtual RealGradient precomputeQpResidual() override;
  virtual RealGradient precomputeQpJacobian() override;

  RankTwoTensor _I;
};
