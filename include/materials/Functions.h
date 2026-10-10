#pragma once

#include "Material.h"
#include "RankTwoTensor.h"

class Functions : public Material
{
public:
  static InputParameters validParams();
  Functions(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const VariableValue & _theta;
  const VariableGradient & _grad_theta;
  const VariableValue & _thetaP;
  const VariableValue & _thetaM;

  MaterialProperty<Real> & _B;
  MaterialProperty<Real> & _R;
  MaterialProperty<Real> & _w;
  MaterialProperty<Real> & _c;
  MaterialProperty<RealVectorValue> & _dB;
  MaterialProperty<RealVectorValue> & _dw;
  MaterialProperty<Real> & _dc;
  MaterialProperty<Real> & _psi;
  MaterialProperty<Real> & _M;
  MaterialProperty<Real> & _energy;
  MaterialProperty<RankTwoTensor> & _d2B;
  MaterialProperty<RankTwoTensor> & _d2w;
  MaterialProperty<Real> & _d2c;
  
  const bool _const_M;
  const Real _M0;
  const Real _alpha;
  const Real _beta;
  const Real _u0;
  const Real _u1_0;
  const Real _B0;
  const int _n;
  const int _m;
  const Real _eps_m;
  const Real _K;

  RankTwoTensor _U;
  RankTwoTensor _I;
};
