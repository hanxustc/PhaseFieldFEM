#pragma once
#include "Kernel.h"
#include "RankTwoTensor.h"

class Bulk : public Kernel
{
public:
  static InputParameters validParams();
  Bulk(const InputParameters & params);

protected:
  const bool _is_sharp;
  const Real _beta;
  const Real _u0;

  const VariableValue & _theta;
  const VariableGradient & _grad_theta;
  const MaterialProperty<Real> & _B;
  const MaterialProperty<Real> & _R;
  const MaterialProperty<Real> & _w;
  const MaterialProperty<Real> & _c;
  const MaterialProperty<Real> & _dc;
  const MaterialProperty<Real> & _M;
  const MaterialProperty<RealVectorValue> & _dw;
  const MaterialProperty<Real> & _d2c;

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
};
