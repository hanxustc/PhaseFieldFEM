#include "qpNeighbors.h"

registerMooseObject("PhaseFieldApp", qpNeighbors);

InputParameters
qpNeighbors::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addRequiredCoupledVar("theta", "The global orientation variable");
  params.addRequiredParam<UserObjectName>("initial_theta", "initialTheta nodal object");
  params.addRequiredParam<Real>("u0", "Lower bound for GBs");
  params.addRequiredParam<int>("every", "Timestep interval for updating P/M neighbors");
  params.addRequiredParam<Real>("dmax", "Maximal sampling distance");
  params.addRequiredParam<Real>("dmin", "Minimal sampling distance");
  params.addParam<bool>("use_smooth", "Use smoothed gradient for neighbor search direction");

  return params;
}

qpNeighbors::qpNeighbors(const InputParameters & params)
  : ElementUserObject(params),
    _theta(adCoupledValue("theta")),
    _grad_theta(adCoupledGradient("theta")),
    _initial_theta(getUserObject<initialTheta>("initial_theta")),
    _theta_init(_initial_theta.getThetaInit()),
    _every(getParam<int>("every")),
    _u0(getParam<Real>("u0")),
    _dmax(getParam<Real>("dmax")),
    _dmin(getParam<Real>("dmin")),
    _use_smooth(getParam<bool>("use_smooth"))
{
  const auto & mesh = _mesh.getMesh();
  libMesh::BoundingBox bbox = libMesh::MeshTools::create_bounding_box(mesh);

  _Lx = bbox.max()(0) - bbox.min()(0);
  _Ly = bbox.max()(1) - bbox.min()(1);

  _h = std::numeric_limits<Real>::max();
  for (const auto & elem : mesh.active_element_ptr_range())
    _h = std::min(_h, elem->hmin());
  _communicator.min(_h);

  if (std::fmod(_dmin, _h) < _dmin / 2.0)
    n = static_cast<int>(_dmin / _h);
  else
    n = static_cast<int>(_dmin / _h) + 1;
}

void
qpNeighbors::meshChanged()
{
  _mesh_changed = true;
  _all_qps.clear();
}

Real
qpNeighbors::periodDist(const libMesh::Point & p1, const libMesh::Point & p2) const
{
  Real dx = std::abs(p1(0) - p2(0));
  Real dy = std::abs(p1(1) - p2(1));
  dx = std::min(dx, _Lx-dx);
  dy = std::min(dy, _Ly-dy);

  return dx*dx + dy*dy;
}

Real
qpNeighbors::Dist(const libMesh::Point & p1, const libMesh::Point & p2) const
{
  Real dx = std::abs(p1(0) - p2(0));
  Real dy = std::abs(p1(1) - p2(1));

  return dx*dx + dy*dy;
}

libMesh::VectorValue<Real>
qpNeighbors::periodDirect(const libMesh::Point & p0, const libMesh::Point & p1) const
{
  Real dx = p1(0) - p0(0);
  Real dy = p1(1) - p0(1);
  if (dx >  0.5*_Lx) dx -= _Lx;
  if (dx < -0.5*_Lx) dx += _Lx;
  if (dy >  0.5*_Ly) dy -= _Ly;
  if (dy < -0.5*_Ly) dy += _Ly;

  return libMesh::VectorValue<Real>(dx, dy, 0.0);
}

libMesh::VectorValue<Real>
qpNeighbors::Direct(const libMesh::Point & p0, const libMesh::Point & p1) const
{
  Real dx = p1(0) - p0(0);
  Real dy = p1(1) - p0(1);

  return libMesh::VectorValue<Real>(dx, dy, 0.0);
}

void
qpNeighbors::buildGrid(Real cutoff)
{
  _grid_x0 = 0.0;
  _grid_y0 = 0.0;
  for (const auto & q : _all_qps) {
    if (q.point(0) < _grid_x0) _grid_x0 = q.point(0);
    if (q.point(1) < _grid_y0) _grid_y0 = q.point(1);
  }

  _grid_dx = cutoff;
  _grid_dy = cutoff;
  _grid_nx = std::max(1, static_cast<int>(std::ceil(_Lx / _grid_dx)));
  _grid_ny = std::max(1, static_cast<int>(std::ceil(_Ly / _grid_dy)));
  _grid_dx = _Lx / _grid_nx;
  _grid_dy = _Ly / _grid_ny;

  _grid.assign(static_cast<std::size_t>(_grid_nx) * _grid_ny, std::vector<std::size_t>());
  for (std::size_t i = 0; i < _all_qps.size(); i++) {
    int cx = static_cast<int>(((_all_qps[i].point(0) - _grid_x0) / _grid_dx));
    int cy = static_cast<int>(((_all_qps[i].point(1) - _grid_y0) / _grid_dy));
    if (cx >= _grid_nx) cx = _grid_nx - 1;
    if (cy >= _grid_ny) cy = _grid_ny - 1;
    if (cx < 0) cx = 0;
    if (cy < 0) cy = 0;
    _grid[static_cast<std::size_t>(cy)*_grid_nx + cx].push_back(i);
  }
}

void
qpNeighbors::initialize()
{
  if (_theta_init.size() > 1)
    _all_qps.clear();
}

void
qpNeighbors::execute()
{
  if (_theta_init.size() > 1) {
    for (std::size_t qp = 0; qp < _qrule->n_points(); qp++) {
      QpInfo info;
      info.elem_id = _current_elem->id();
      info.qp_id = qp;
      info.point = _q_point[qp];
      info.theta = MetaPhysicL::raw_value(_theta[qp]);
      info.grad_theta = MetaPhysicL::raw_value(_grad_theta[qp]);
      _all_qps.push_back(info);
    }
  }
}

void
qpNeighbors::finalize()
{
  if (_theta_init.size() > 1) {
    if (!_mesh.getMesh().is_serial()) {
      const std::size_t stride = 9;
      std::vector<Real> buf;
      buf.reserve(_all_qps.size() * stride);
      for (const auto & q : _all_qps) {
        buf.push_back(static_cast<Real>(q.elem_id));
        buf.push_back(static_cast<Real>(q.qp_id));
        buf.push_back(q.point(0));
        buf.push_back(q.point(1));
        buf.push_back(q.point(2));
        buf.push_back(q.theta);
        buf.push_back(q.grad_theta(0));
        buf.push_back(q.grad_theta(1));
        buf.push_back(q.grad_theta(2));
      }
      _communicator.allgather(buf);

      _all_qps.clear();
      _all_qps.reserve(buf.size() / stride);
      for (std::size_t i = 0; i < buf.size(); i += stride) {
        QpInfo q;
        q.elem_id = static_cast<unsigned int>(buf[i+0]);
        q.qp_id = static_cast<unsigned int>(buf[i+1]);
        q.point = libMesh::Point(buf[i+2], buf[i+3], buf[i+4]);
        q.theta = buf[i+5];
        q.grad_theta = libMesh::VectorValue<Real>(buf[i+6], buf[i+7], buf[i+8]);
        _all_qps.push_back(q);
      }
    }

    const Real eps = 1e-6;
    _rebuilt = false;

    if (_t_step % _every == 0 || _mesh_changed) {
      const Real smooth_cutoff = std::sqrt(6.0) * _dmin;
      const Real grid_cutoff = std::max(smooth_cutoff, _dmax);
      buildGrid(grid_cutoff);
      const int search_r_smooth = std::max(1, static_cast<int>(std::ceil(smooth_cutoff / _grid_dx)));
      const int search_r_neigh  = std::max(1, static_cast<int>(std::ceil(_dmax / _grid_dx)));

      if (_use_smooth) {
        const Real sigma2 = _dmin * _dmin;
        const Real cutoff2 = 6.0 * sigma2;
        _smooth_grad.clear();
        _smooth_grad.reserve(_all_qps.size());

        for (const auto & q0 : _all_qps) {
          int cx0 = static_cast<int>((q0.point(0) - _grid_x0) / _grid_dx);
          int cy0 = static_cast<int>((q0.point(1) - _grid_y0) / _grid_dy);
          if (cx0 >= _grid_nx) cx0 = _grid_nx - 1;
          if (cy0 >= _grid_ny) cy0 = _grid_ny - 1;
          if (cx0 < 0) cx0 = 0;
          if (cy0 < 0) cy0 = 0;

          Real wsum = 0.0, gx = 0.0, gy = 0.0;
          for (int dy = -search_r_smooth; dy <= search_r_smooth; dy++)
            for (int dx = -search_r_smooth; dx <= search_r_smooth; dx++) {
              int cx1 = (cx0 + dx % _grid_nx + _grid_nx) % _grid_nx;
              int cy1 = (cy0 + dy % _grid_ny + _grid_ny) % _grid_ny;
              for (const auto & i1 : _grid[static_cast<std::size_t>(cy1) * _grid_nx + cx1]) {
                const auto & q1 = _all_qps[i1];
                const Real d2 = periodDist(q0.point, q1.point);
                if (d2 > cutoff2) continue;
                const Real w = std::exp(-d2 / sigma2);
                gx += w * q1.grad_theta(0);
                gy += w * q1.grad_theta(1);
                wsum += w;
              }
            }
          _smooth_grad[key(q0.elem_id, q0.qp_id)] = (wsum > 1e-16) ? libMesh::VectorValue<Real>(gx/wsum, gy/wsum, 0.0) : q0.grad_theta;
        }
      }
      _neigh_plus_qps.clear();
      _neigh_plus_qps.reserve(_all_qps.size());
      _neigh_minus_qps.clear();
      _neigh_minus_qps.reserve(_all_qps.size());

      for (std::size_t i0 = 0; i0 < _all_qps.size(); i0++) {
        const auto & q0 = _all_qps[i0];
        Real ref_norm;
        libMesh::VectorValue<Real> ref_grad;
        if (_use_smooth) {
          const auto & grad_smooth = _smooth_grad[key(q0.elem_id, q0.qp_id)];
          ref_norm = std::sqrt(grad_smooth(0)*grad_smooth(0) + grad_smooth(1)*grad_smooth(1));
          ref_grad = grad_smooth;
        } else {
          ref_norm = q0.grad_theta.norm();
          ref_grad = q0.grad_theta;
        }

        if (ref_norm <= _u0)
          continue;

        auto & neighP = _neigh_plus_qps[key(q0.elem_id, q0.qp_id)];
        auto & neighM = _neigh_minus_qps[key(q0.elem_id, q0.qp_id)];
        neighP.clear();
        neighM.clear();

        libMesh::VectorValue<Real> grad;
        grad = ref_grad / ref_norm;
        std::vector<std::pair<Real, std::size_t>> plus_list;
        std::vector<std::pair<Real, std::size_t>> minus_list;

        int cx0 = static_cast<int>((q0.point(0) - _grid_x0) / _grid_dx);
        int cy0 = static_cast<int>((q0.point(1) - _grid_y0) / _grid_dy);
        if (cx0 >= _grid_nx) cx0 = _grid_nx - 1;
        if (cy0 >= _grid_ny) cy0 = _grid_ny - 1;
        if (cx0 < 0) cx0 = 0;
        if (cy0 < 0) cy0 = 0;

        Real threshold = 0.999;
        while (true) {
          plus_list.clear();
          minus_list.clear();

          for (int dy = -search_r_neigh; dy <= search_r_neigh; dy++)
            for (int dx = -search_r_neigh; dx <= search_r_neigh; dx++) {
              int cx1 = (cx0 + dx % _grid_nx + _grid_nx) % _grid_nx;
              int cy1 = (cy0 + dy % _grid_ny + _grid_ny) % _grid_ny;
              for (const auto & i1 : _grid[static_cast<std::size_t>(cy1) * _grid_nx + cx1]) {
                if (i1 == i0)
                  continue;
                const auto & q1 = _all_qps[i1];

                // libMesh::VectorValue<Real> v = Direct(q0.point, q1.point);
                libMesh::VectorValue<Real> v = periodDirect(q0.point, q1.point);
                const Real r2 = v(0)*v(0) + v(1)*v(1);
                if (r2 > _dmax*_dmax || r2 < 0.5*_dmax*_dmax)
                  continue;

                const Real r = std::sqrt(r2);
                const Real proj = v(0)*grad(0) + v(1)*grad(1);
                const Real cosine = proj / r;

                if (std::abs(cosine) < threshold)
                  continue;

                if (proj > 0)
                  plus_list.emplace_back(r2, i1);
                else if (proj < 0)
                  minus_list.emplace_back(r2, i1);
              }
            }

          if (!plus_list.empty() && !minus_list.empty())
            break;
          if (threshold <= 0.0)
            break;
          threshold -= 0.003;
        }

        auto r2_comp = [](const std::pair<Real, std::size_t> &a,
                          const std::pair<Real, std::size_t> &b)
        { return a.first < b.first; };

        std::sort(plus_list.begin(), plus_list.end(), r2_comp);
        std::sort(minus_list.begin(), minus_list.end(), r2_comp);

        if (plus_list.empty() && minus_list.empty()) {
          neighP.push_back(i0);
          neighM.push_back(i0);
          continue;
        }

        if (plus_list.empty()) {
          for (std::size_t k = 0; k < minus_list.size(); k++)
            plus_list.emplace_back(0.0, i0);
        }

        if (minus_list.empty()) {
          for (std::size_t k = 0; k < plus_list.size(); k++)
            minus_list.emplace_back(0.0, i0);
        }

        const std::size_t n_layers = std::max(plus_list.size(), minus_list.size());
        neighP.reserve(n_layers);
        neighM.reserve(n_layers);
        for (std::size_t layer = 0; layer < n_layers; layer++) {
          neighP.push_back(plus_list[std::min(layer, plus_list.size()-1)].second);
          neighM.push_back(minus_list[std::min(layer, minus_list.size()-1)].second);
        }
      }

      _mesh_changed = false;
      _rebuilt = true;

      if (processor_id() == 0)
        std::cout << "\nP/M neighbors rebuilt at t_step = "
                  << _t_step << std::endl;
    }
  }
}

void
qpNeighbors::threadJoin(const UserObject & other)
{
  if (_theta_init.size() > 1) {
    const qpNeighbors & rhs = dynamic_cast<const qpNeighbors &>(other);
    _all_qps.insert(_all_qps.end(), rhs._all_qps.begin(), rhs._all_qps.end());
  }
}

const std::vector<qpNeighbors::QpInfo> &
qpNeighbors::getAllQps() const
{
  return _all_qps;
}

const std::vector<std::size_t> &
qpNeighbors::getPlusNeighbors(dof_id_type elem_id, unsigned int qp_id) const
{
  static const std::vector<std::size_t> empty;
  const auto it = _neigh_plus_qps.find(key(elem_id, qp_id));
  if (it != _neigh_plus_qps.end())
    return it->second;
  return empty;
}

const std::vector<std::size_t> &
qpNeighbors::getMinusNeighbors(dof_id_type elem_id, unsigned int qp_id) const
{
  static const std::vector<std::size_t> empty;
  const auto it = _neigh_minus_qps.find(key(elem_id, qp_id));
  if (it != _neigh_minus_qps.end())
    return it->second;
  return empty;
}

Real
qpNeighbors::getSmoothGradNorm(dof_id_type elem_id, unsigned int qp_id) const
{
  const auto it = _smooth_grad.find(key(elem_id, qp_id));
  if (it != _smooth_grad.end()) {
    const auto & g = it->second;
    return std::sqrt(g(0)*g(0) + g(1)*g(1));
  }
  return 0.0;
}

const libMesh::VectorValue<Real> &
qpNeighbors::getSmoothGrad(dof_id_type elem_id, unsigned int qp_id) const
{
  static const libMesh::VectorValue<Real> zero(0.0, 0.0, 0.0);
  const auto it = _smooth_grad.find(key(elem_id, qp_id));
  if (it != _smooth_grad.end())
    return it->second;
  return zero;
}