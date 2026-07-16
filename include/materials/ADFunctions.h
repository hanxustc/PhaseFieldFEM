#pragma once

#include "Material.h"
#include "ADMaterial.h"

class ADFunctions : public ADMaterial
{
public:
  static InputParameters validParams();
  ADFunctions(const InputParameters & parameters);
  /*
  struct OtherQpInfo
  {
    unsigned int elem_id;
    unsigned int qp_id;
    libMesh::Point point;
    Real theta;
  };
*/
protected:
  virtual void computeQpProperties() override;

  // Coupled variable
  const ADVariableValue & _theta;
  const ADVariableGradient & _grad_theta;
  const VariableValue & _thetaP;
  const VariableValue & _thetaM;
  
  // Material functions
  ADMaterialProperty<Real> & _B;
  ADMaterialProperty<Real> & _R;
  ADMaterialProperty<Real> & _w;
  ADMaterialProperty<Real> & _c;
  ADMaterialProperty<RealVectorValue> & _dB;
  ADMaterialProperty<RealVectorValue> & _dw;
  ADMaterialProperty<Real> & _dc;
  ADMaterialProperty<Real> & _psi;
  ADMaterialProperty<Real> & _M;
  ADMaterialProperty<Real> & _energy;
  
  // Input parameters
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
};