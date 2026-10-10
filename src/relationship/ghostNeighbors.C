#include "ghostNeighbors.h"
#include "MooseMesh.h"
#include "MooseApp.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include <cmath>
#include <unordered_set>

registerMooseObject("PhaseFieldApp", ghostNeighbors);

using namespace libMesh;

InputParameters
ghostNeighbors::validParams()
{
  InputParameters params = RelationshipManager::validParams();
  params.addRequiredParam<Real>("distance", "Physical radius to ghost out");
  params.set<bool>("attach_geometric_early") = false;
  return params;
}

ghostNeighbors::ghostNeighbors(const InputParameters & params)
  : RelationshipManager(params), _distance(getParam<Real>("distance")) {}

ghostNeighbors::ghostNeighbors(const ghostNeighbors & other)
  : RelationshipManager(other), _distance(other._distance), _Lx(other._Lx), _Ly(other._Ly), _xmin(other._xmin), _ymin(other._ymin) {}

void
ghostNeighbors::internalInitWithMesh(const MeshBase & mesh)
{
  const libMesh::BoundingBox bbox = MeshTools::create_bounding_box(mesh);
  _xmin = bbox.min()(0);
  _ymin = bbox.min()(1);
  _Lx = bbox.max()(0) - _xmin;
  _Ly = bbox.max()(1) - _ymin;
}

Real
ghostNeighbors::periodDist2(const libMesh::Point & a, const libMesh::Point & b) const
{
  Real dx = std::abs(a(0) - b(0));
  Real dy = std::abs(a(1) - b(1));

  if (_Lx > 0.0)
    dx = std::min(dx, _Lx - dx);
  if (_Ly > 0.0)
    dy = std::min(dy, _Ly - dy);
  return dx * dx + dy * dy;
}

void
ghostNeighbors::operator()(const MeshBase::const_element_iterator & range_begin, const MeshBase::const_element_iterator & range_end, const processor_id_type p, map_type & coupled_elements)
{
  if (range_begin == range_end)
    return;

  mooseAssert(_mesh, "Mesh must be set");

  const dof_id_type nelem_now = _mesh->n_elem();
  if (_hmax_all < 0.0 || _hmax_nelem != nelem_now) {
    Real h = 0.0;
    for (const Elem * const elem : _mesh->active_element_ptr_range())
      h = std::max(h, elem->hmax());
    _hmax_all = h;
    _hmax_nelem = nelem_now;
  }
  Real hmax = _hmax_all;
  for (auto it = range_begin; it != range_end; ++it)
    hmax = std::max(hmax, (*it)->hmax());

  const Real reach = _distance + hmax;
  const Real target = reach / 6.0;
  const int nx = std::max(1, static_cast<int>(std::ceil(_Lx / target)));
  const int ny = std::max(1, static_cast<int>(std::ceil(_Ly / target)));
  const Real bx = _Lx / nx;
  const Real by = _Ly / ny;

  auto bin_of = [&](const libMesh::Point & pt, int & cx, int & cy)
  {
    cx = static_cast<int>((pt(0) - _xmin) / bx);
    cy = static_cast<int>((pt(1) - _ymin) / by);
    cx = std::min(std::max(cx, 0), nx - 1);
    cy = std::min(std::max(cy, 0), ny - 1);
  };

  const std::size_t nbin = static_cast<std::size_t>(nx) * ny;
  std::vector<char> occupied(nbin, 0);
  std::vector<int> occ_list;
  std::vector<std::vector<libMesh::Point>> pts(nbin);
  for (auto it = range_begin; it != range_end; ++it) {
    const libMesh::Point ctr = (*it)->vertex_average();
    int cx, cy;
    bin_of(ctr, cx, cy);
    const std::size_t b = static_cast<std::size_t>(cy) * nx + cx;
    if (!occupied[b]) {
      occupied[b] = 1;
      occ_list.push_back(static_cast<int>(b));
    }
    pts[b].push_back(ctr);
  }

  const Real dilate = reach + std::sqrt(bx*bx + by*by);
  const int rx = static_cast<int>(std::ceil(dilate / bx)) + 1;
  const int ry = static_cast<int>(std::ceil(dilate / by)) + 1;
  std::vector<char> mask(nbin, 0);
  for (const int b : occ_list) {
    const int cx = b % nx, cy = b / nx;
    for (int dy = -ry; dy <= ry; dy++)
      for (int dx = -rx; dx <= rx; dx++) {
        const Real ddx = dx * bx, ddy = dy * by;
        if (ddx*ddx + ddy*ddy > dilate*dilate) continue;
        const int qx = ((cx + dx) % nx + nx) % nx;
        const int qy = ((cy + dy) % ny + ny) % ny;
        mask[static_cast<std::size_t>(qy) * nx + qx] = 1;
      }
  }

  const int sx = static_cast<int>(std::ceil(reach / bx)) + 1;
  const int sy = static_cast<int>(std::ceil(reach / by)) + 1;
  const Real reach2 = reach * reach;
  static const CouplingMatrix * const null_mat = nullptr;
  for (const Elem * const elem : _mesh->active_element_ptr_range()) {
    if (elem->processor_id() == p)
      continue;
    const libMesh::Point ctr = elem->vertex_average();
    int cx, cy;
    bin_of(ctr, cx, cy);
    if (!mask[static_cast<std::size_t>(cy) * nx + cx])
      continue;
    bool within = false;
    for (int dy = -sy; dy <= sy && !within; dy++)
      for (int dx = -sx; dx <= sx && !within; dx++) {
        const int qx = ((cx + dx) % nx + nx) % nx;
        const int qy = ((cy + dy) % ny + ny) % ny;
        for (const auto & q : pts[static_cast<std::size_t>(qy) * nx + qx])
          if (periodDist2(ctr, q) <= reach2) {
            within = true;
            break;
          }
      }
    if (within) {
      coupled_elements.emplace(elem, null_mat);
    }
  }
}

std::unique_ptr<GhostingFunctor>
ghostNeighbors::clone() const
{
  return _app.getFactory().copyConstruct(*this);
}

std::string
ghostNeighbors::getInfo() const
{
  std::ostringstream oss;
  oss << "ghostNeighbors (radius " << _distance << ')';
  return oss.str();
}

bool
ghostNeighbors::operator>=(const RelationshipManager & other) const
{
  const auto * rm = dynamic_cast<const ghostNeighbors *>(&other);
  if (!rm)
    return false;
  return _distance >= rm->_distance && baseGreaterEqual(*rm);
}
