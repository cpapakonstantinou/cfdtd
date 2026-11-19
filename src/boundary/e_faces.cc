/** \file   boundary/e_faces.cc
 *  \brief  psi updates of the E faces of the CPML, inc/boundary/cpml.h.
 */
#include "boundary/cpml.h"
#include "materials/models.h"

#include "materials/stretch.h"
#include <algorithm>
#include "grid/lattice.h"
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::fields::Hx;
using fdtd::fields::Hy;
using fdtd::fields::Hz;
using fdtd::fields::ex_nl;
using fdtd::fields::ey_nl;
using fdtd::fields::ez_nl;
using fdtd::fields::hx_nl;
using fdtd::fields::hy_nl;
using fdtd::fields::hz_nl;
using fdtd::grid::subdomain;
using fdtd::materials::max_slots;
using fdtd::materials::matind_x_full;
using fdtd::materials::matind_y_full;
using fdtd::materials::matind_z_full;
using fdtd::materials::background_y;

namespace
{
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace

namespace
{
fdtd::boundary::Cpml &the_cpml = fdtd::boundary::cpml(); ///< cached: the accessor is out of line
} // namespace

namespace fdtd::boundary
{
namespace
{
enum class Axis { X, Y, Z };

/**
 * @brief One CPML face of one E component, and of its separate nonlinear field
 * @verbatim
      rows (i, j) in [i0, i1) x [j0, j1), cells k in [k0, k1)
      grading index s:  near face  n + pad_nea      far face  cells - n + pad_far    (n = i, j or k)
      far faces store psi from the face start, i, j or k minus *_far_cpml_start, and skip index 0
      factor: sign * ce_2_d of a Drude or Lorentz background; other dispersive: none,
              sign * ce_2 of a dielectric background, or of the cell's material with mode_source
      background: first / second medium for the y faces, the y layer's medium otherwise
   @endverbatim
 * @param d  d(h, i, j, k): the derivative of the H array h for this face
 */
template <class D>
void e_face(Axis axis, bool far, int i0, int i1, int j0, int j1, int k0, int k1, double sign,
		   view3<double> e, view3<double> e_nl, view3<int> full, view3<double> psi, view3<double> psi_nl,
		   const array1<double> &c1, const array1<double> &c2, view3<double> h, view3<double> h_nl, D &&d)
{
	const auto &mp = fdtd::materials::model_parameters();
	const int cells = axis == Axis::X ? subdomain.cells_x : axis == Axis::Y ? subdomain.cells_y : subdomain.cells_z;
	const int pad = axis == Axis::X ? (far ? subdomain.cpml_pad_far_x : subdomain.cpml_pad_nea_x)
							  : axis == Axis::Y ? (far ? subdomain.cpml_pad_far_y : subdomain.cpml_pad_nea_y)
												: (far ? subdomain.cpml_pad_far_z : subdomain.cpml_pad_nea_z);
	const int extra = axis == Axis::X ? subdomain.i_far_cpml_start : axis == Axis::Y ? subdomain.j_far_cpml_start : subdomain.k_far_cpml_start;
	const int ka = (far && axis == Axis::Z) ? std::max(k0, 1) : k0;
	const int off = (far && axis == Axis::Z) ? extra : 0;
	for (int i = i0; i < i1; ++i)
	{
		if (far && axis == Axis::X && i == 0)
			continue;
		for (int j = j0; j < j1; ++j)
		{
			if (far && axis == Axis::Y && j == 0)
				continue;
			const int m = axis == Axis::Y ? (far ? the_cpml.background.far_y : the_cpml.background.near_y) : background_y[j + subdomain.j_for_uni];
			const int pi = (far && axis == Axis::X) ? i - extra : i;
			const int pj = (far && axis == Axis::Y) ? j - extra : j;
			auto s = [&](int k) { const int n = axis == Axis::X ? i : axis == Axis::Y ? j : k; return far ? cells - n + pad : n + pad; };
			auto row = [&](view3<double> f, view3<double> p, view3<double> hh, auto &&factor) {
				fdtd::materials::stretch::psi_update(f[i][j], p[pi][pj], off, ka, k1, c1, c2, s, factor, [&](int k) { return d(hh, i, j, k); });
			};
			if (fdtd::materials::dispersive(m))
			{
				if (fdtd::materials::adds_terms_after(m))
				{
					const double factor = sign * fdtd::materials::after_terms_ce2(m);
					row(e, psi, h, [=](int) { return factor; });
				}
			}
			else if (!the_cpml.background.cellwise)
			{
				const double factor = sign * mp.ce_2[m];
				row(e, psi, h, [=](int) { return factor; });
				if (the_cpml.harmonic)
					row(e_nl, psi_nl, h_nl, [=](int) { return factor; });
			}
			else
				row(e, psi, h, [&](int k) { return sign * mp.ce_2[full[i][j][k]]; });
		}
	}
}

/// end of a far face: the last cell when the far CPML is padded on this rank, else one before
int far_end(int cells, int pad_far) { return pad_far > 0 ? cells + 1 : cells; }
} // namespace

/// CPML auxiliary fields of Ex: y faces (+), z faces (-)
void update_psi_ex()
{
	const auto &S = subdomain;
	auto dHz_y = [](view3<double> h, int i, int j, int k) { return (h[i][j][k] - h[i][j - 1][k]) * lat.inv_d[1]; };
	auto dHy_z = [](view3<double> h, int i, int j, int k) { return (h[i][j][k] - h[i][j][k - 1]) * lat.inv_d[2]; };
	const int nz = S.owned_cells_z + 1;
	if (subdomain.nea_cpml_y())
		e_face(Axis::Y, false, 0, S.cells_x, 1, S.j_nea_cpml_end + 1, 1, nz, +1, Ex, ex_nl, matind_x_full, the_cpml.psi[0][0][0][1][0], the_cpml.psi[1][0][0][1][0], the_cpml.profile[1][0].b[0], the_cpml.profile[1][0].c[0], Hz, hz_nl, dHz_y);
	if (subdomain.far_cpml_y())
		e_face(Axis::Y, true, 0, S.cells_x, S.j_far_cpml_start, far_end(S.cells_y, S.cpml_pad_far_y), 1, nz, +1, Ex, ex_nl, matind_x_full, the_cpml.psi[0][0][0][1][1], the_cpml.psi[1][0][0][1][1], the_cpml.profile[1][1].b[0], the_cpml.profile[1][1].c[0], Hz, hz_nl, dHz_y);
	if (!lat.periodic_z && subdomain.nea_cpml_z())
		e_face(Axis::Z, false, 0, S.cells_x, 1, S.owned_cells_y + 1, 1, S.k_nea_cpml_end + 1, -1, Ex, ex_nl, matind_x_full, the_cpml.psi[0][0][0][2][0], the_cpml.psi[1][0][0][2][0], the_cpml.profile[2][0].b[0], the_cpml.profile[2][0].c[0], Hy, hy_nl, dHy_z);
	if (!lat.periodic_z && subdomain.far_cpml_z())
		e_face(Axis::Z, true, 0, S.cells_x, 1, S.owned_cells_y + 1, S.k_far_cpml_start, far_end(S.cells_z, S.cpml_pad_far_z), -1, Ex, ex_nl, matind_x_full, the_cpml.psi[0][0][0][2][1], the_cpml.psi[1][0][0][2][1], the_cpml.profile[2][0].b[0], the_cpml.profile[2][0].c[0], Hy, hy_nl, dHy_z);
}

/// CPML auxiliary fields of Ey: x faces (-), z faces (+)
void update_psi_ey()
{
	const auto &S = subdomain;
	auto dHz_x = [](view3<double> h, int i, int j, int k) { return (h[i][j][k] - h[i - 1][j][k]) * lat.inv_d[0]; };
	auto dHx_z = [](view3<double> h, int i, int j, int k) { return (h[i][j][k] - h[i][j][k - 1]) * lat.inv_d[2]; };
	const int nz = S.owned_cells_z + 1;
	if (!lat.periodic_x && subdomain.nea_cpml_x())
		e_face(Axis::X, false, 1, S.i_nea_cpml_end + 1, 0, S.cells_y, 1, nz, -1, Ey, ey_nl, matind_y_full, the_cpml.psi[0][0][1][0][0], the_cpml.psi[1][0][1][0][0], the_cpml.profile[0][0].b[0], the_cpml.profile[0][0].c[0], Hz, hz_nl, dHz_x);
	if (!lat.periodic_x && subdomain.far_cpml_x())
		e_face(Axis::X, true, S.i_far_cpml_start, far_end(S.cells_x, S.cpml_pad_far_x), 0, S.cells_y, 1, nz, -1, Ey, ey_nl, matind_y_full, the_cpml.psi[0][0][1][0][1], the_cpml.psi[1][0][1][0][1], the_cpml.profile[0][0].b[0], the_cpml.profile[0][0].c[0], Hz, hz_nl, dHz_x);
	if (!lat.periodic_z && subdomain.nea_cpml_z())
		e_face(Axis::Z, false, 1, S.owned_cells_x + 1, 0, S.cells_y, 1, S.k_nea_cpml_end + 1, +1, Ey, ey_nl, matind_y_full, the_cpml.psi[0][0][1][2][0], the_cpml.psi[1][0][1][2][0], the_cpml.profile[2][0].b[0], the_cpml.profile[2][0].c[0], Hx, hx_nl, dHx_z);
	if (!lat.periodic_z && subdomain.far_cpml_z())
		e_face(Axis::Z, true, 1, S.owned_cells_x + 1, 0, S.cells_y, S.k_far_cpml_start, far_end(S.cells_z, S.cpml_pad_far_z), +1, Ey, ey_nl, matind_y_full, the_cpml.psi[0][0][1][2][1], the_cpml.psi[1][0][1][2][1], the_cpml.profile[2][0].b[0], the_cpml.profile[2][0].c[0], Hx, hx_nl, dHx_z);
}

/// CPML auxiliary fields of Ez: x faces (+), y faces (-)
void update_psi_ez()
{
	const auto &S = subdomain;
	auto dHy_x = [](view3<double> h, int i, int j, int k) { return (h[i][j][k] - h[i - 1][j][k]) * lat.inv_d[0]; };
	auto dHx_y = [](view3<double> h, int i, int j, int k) { return (h[i][j][k] - h[i][j - 1][k]) * lat.inv_d[1]; };
	const int nz = S.cells_z;
	if (!lat.periodic_x && subdomain.nea_cpml_x())
		e_face(Axis::X, false, 1, S.i_nea_cpml_end + 1, 1, S.owned_cells_y + 1, 0, nz, +1, Ez, ez_nl, matind_z_full, the_cpml.psi[0][0][2][0][0], the_cpml.psi[1][0][2][0][0], the_cpml.profile[0][0].b[0], the_cpml.profile[0][0].c[0], Hy, hy_nl, dHy_x);
	if (!lat.periodic_x && subdomain.far_cpml_x())
		e_face(Axis::X, true, S.i_far_cpml_start, far_end(S.cells_x, S.cpml_pad_far_x), 1, S.owned_cells_y + 1, 0, nz, +1, Ez, ez_nl, matind_z_full, the_cpml.psi[0][0][2][0][1], the_cpml.psi[1][0][2][0][1], the_cpml.profile[0][0].b[0], the_cpml.profile[0][0].c[0], Hy, hy_nl, dHy_x);
	if (subdomain.nea_cpml_y())
		e_face(Axis::Y, false, 1, S.owned_cells_x + 1, 1, S.j_nea_cpml_end + 1, 0, nz, -1, Ez, ez_nl, matind_z_full, the_cpml.psi[0][0][2][1][0], the_cpml.psi[1][0][2][1][0], the_cpml.profile[1][0].b[0], the_cpml.profile[1][0].c[0], Hx, hx_nl, dHx_y);
	if (subdomain.far_cpml_y())
		e_face(Axis::Y, true, 1, S.owned_cells_x + 1, S.j_far_cpml_start, far_end(S.cells_y, S.cpml_pad_far_y), 0, nz, -1, Ez, ez_nl, matind_z_full, the_cpml.psi[0][0][2][1][1], the_cpml.psi[1][0][2][1][1], the_cpml.profile[1][1].b[0], the_cpml.profile[1][1].c[0], Hx, hx_nl, dHx_y);
}
} // namespace fdtd::boundary
