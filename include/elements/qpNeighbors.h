#pragma once

#include "initialTheta.h"
#include "ElementUserObject.h"
#include "libmesh/point.h"
#include "libmesh/utility.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_tools.h"
#include <unordered_map>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>

class qpNeighbors : public ElementUserObject
{
public:
  struct QpInfo
  {
    unsigned int elem_id;
    unsigned int qp_id;
    libMesh::Point point;
    Real theta;
    libMesh::VectorValue<Real> grad_theta;
  };

  static InputParameters validParams();
  qpNeighbors(const InputParameters &parameters);
  
  std::unordered_map<std::uint64_t, std::vector<std::size_t>> _neigh_plus_qps;
  std::unordered_map<std::uint64_t, std::vector<std::size_t>> _neigh_minus_qps;
  std::unordered_map<std::uint64_t, libMesh::VectorValue<Real>> _smooth_grad;

  std::vector<QpInfo> _all_qps;
  const std::vector<std::size_t> &getPlusNeighbors(dof_id_type elem_id, unsigned int qp_id) const;
  const std::vector<std::size_t> &getMinusNeighbors(dof_id_type elem_id, unsigned int qp_id) const;
  const std::vector<QpInfo> &getAllQps() const;

  libMesh::VectorValue<Real> periodDirect(const libMesh::Point &p0, const libMesh::Point &p1) const;
  libMesh::VectorValue<Real> Direct(const libMesh::Point &p0, const libMesh::Point &p1) const;

  bool neighborsRebuilt() const { return _rebuilt; }
  bool useSmooth() const { return _use_smooth; }
  Real getSmoothGradNorm(dof_id_type elem_id, unsigned int qp_id) const;
  const libMesh::VectorValue<Real> & getSmoothGrad(dof_id_type elem_id, unsigned int qp_id) const;

protected:
  virtual void execute() override;
  virtual void initialize() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & other) override;
  virtual void meshChanged() override;

  Real periodDist(const libMesh::Point &p1, const libMesh::Point &p2) const;
  Real Dist(const libMesh::Point &p1, const libMesh::Point &p2) const;
  static inline std::uint64_t key(dof_id_type elem_id, unsigned int qp_id) {
    return (std::uint64_t(elem_id) << 32) | std::uint64_t(qp_id);
  }
  
  const ADVariableValue & _theta;
  const ADVariableGradient & _grad_theta;
  const initialTheta & _initial_theta;
  const std::vector<Real> & _theta_init;
  
  bool _mesh_changed = false;
  bool _rebuilt = false;
  const bool _use_smooth;
  const int _every;
  const Real _u0;
  const Real _dmax;
  const Real _dmin;
  
  Real _Lx = 0.0;
  Real _Ly = 0.0;
  Real _h = 0.0;
  int n = 0;
  
  std::vector<std::vector<std::size_t>> _grid;
  int _grid_nx = 0;
  int _grid_ny = 0;
  Real _grid_dx = 0.0;
  Real _grid_dy = 0.0;
  Real _grid_x0 = 0.0;
  Real _grid_y0 = 0.0;
  void buildGrid(Real cutoff);
};