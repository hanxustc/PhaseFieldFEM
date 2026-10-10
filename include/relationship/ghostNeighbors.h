#pragma once

#include "RelationshipManager.h"
#include "libmesh/bounding_box.h"
#include <vector>

class ghostNeighbors : public RelationshipManager
{
public:
  static InputParameters validParams();

  ghostNeighbors(const InputParameters & params);
  ghostNeighbors(const ghostNeighbors & other);

  void operator()(const MeshBase::const_element_iterator & range_begin,
                  const MeshBase::const_element_iterator & range_end,
                  processor_id_type p,
                  map_type & coupled_elements) override;

  std::unique_ptr<GhostingFunctor> clone() const override;
  std::string getInfo() const override;
  bool operator>=(const RelationshipManager & other) const override;

protected:
  virtual void internalInitWithMesh(const MeshBase &) override;

  const Real _distance;
  Real _Lx = 0.0;
  Real _Ly = 0.0;
  Real _xmin = 0.0;
  Real _ymin = 0.0;

  mutable Real _hmax_all = -1.0;
  mutable dof_id_type _hmax_nelem = 0;
  
  Real periodDist2(const libMesh::Point & a, const libMesh::Point & b) const;
};
