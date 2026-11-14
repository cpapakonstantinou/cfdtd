/**
 * @file sources/tfsf.cc
 * @brief TF/SF corrections of the plane-wave source, sources/tfsf.h
 *
 * One plane loop serves all faces; each face gives the field, the plane, its ranges and the signed
 * correction per cell., Before: twelve functions in planeWave.cc, with dead Bloch and backward-pulse branches
 */
#include "sources/tfsf.h"
#include "materials/models.h"
#include "sources/state.h"
#include "sources/signal.h"

#include "sources/aux_grid.h"
#include "update/state.h"
#include "grid/lattice.h"
#include "sources/spec.h"
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::fields::Hx;
using fdtd::fields::Hy;
using fdtd::fields::Hz;
using fdtd::grid::subdomain;
using fdtd::materials::DispMateOffsetD;
using fdtd::materials::Drude;
using fdtd::materials::DrudeCP;
using fdtd::materials::Lorentz;
using fdtd::materials::NRMAXMATERIALS;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::sources; // transitional: module state, sources/state.h

namespace fdtd::sources
{
namespace
{
/// the 1D source line and the periodic reference line, cached: tfsf_e_terms runs per cell
const AuxLine &src = aux_grid().line[0];
const AuxLine &ref = aux_grid().line[1];

/**
 * @brief f += v(a, b) on the plane n_axis = at
 * @param a0, a1  range [a0, a1) of the first other axis, x for the y and z planes, y for the x plane
 * @param b0, b1  range [b0, b1) of the second other axis, z for the x and y planes, y for the z plane
 */
template <class V>
void correct_plane(view3<double> f, int axis, int at, int a0, int a1, int b0, int b1, V &&v)
{
	for (int a = a0; a < a1; ++a)
		for (int b = b0; b < b1; ++b)
		{
			double &c = axis == 0 ? f[at][a][b] : axis == 1 ? f[a][at][b] : f[a][b][at];
			c += v(a, b);
		}
}

/// E coefficient of the TF/SF correction in background m; false where the model adds it itself
bool background_e_factor(int m, double &c)
{
	const auto &mp = fdtd::materials::model_parameters();
	if (m >= DispMateOffsetD && m < NRMAXMATERIALS)
	{
		const int mm = m - DispMateOffsetD;
		if (mm >= DrudeCP && mm < Lorentz)
		{
			c = (mm < Drude ? fdtd::materials::coeffs<fdtd::materials::models::Drude>(mm).ce2 : fdtd::materials::coeffs<fdtd::materials::models::Lorentz>(mm).ce2);
			return true;
		}
		return false;
	}
	c = mp.ce_2[m];
	return true;
}

bool periodic() { return lat.periodic_x && lat.periodic_z; }
} // namespace

void tfsf_h(Face f)
{
	switch (f)
	{
	case Face::XNea: // Hy from Ez_inc
	case Face::XFar:
	{
		const double s = f == Face::XNea ? -1 : 1;
		correct_plane(Hy, 0, f == Face::XNea ? ipw_near_h : subdomain.ipw_2, subdomain.jpw_min, subdomain.jpw_max + 1, subdomain.kpw_min, subdomain.kpw_max, [&](int j, int) { return s == 1 ? fdtd::sources::chx * src.e[1][j + subdomain.j_for_uni] : -(fdtd::sources::chx * src.e[1][j + subdomain.j_for_uni]); });
		break;
	}
	case Face::ZNea: // Hy from Ex_inc
	case Face::ZFar:
	{
		const bool nea = f == Face::ZNea;
		correct_plane(Hy, 2, nea ? kpw_near_h : subdomain.kpw_2, subdomain.ipw_min, subdomain.ipw_max, subdomain.jpw_min, subdomain.jpw_max + 1, [&](int, int j) { return nea ? fdtd::sources::chz * src.e[0][j + subdomain.j_for_uni] : -(fdtd::sources::chz * src.e[0][j + subdomain.j_for_uni]); });
		break;
	}
	case Face::YNea:
		// with periodic x and z: the reference line; else the 1D grid
		if (zpol)
			correct_plane(Hx, 1, jpw_near_h, subdomain.ipw_min, subdomain.ipw_max + 1, subdomain.kpw_min, subdomain.kpw_max, [&](int /*i*/, int /*k*/) {
				if (periodic())
					return fdtd::sources::chy * ref.e[1][jpw_start_uni];
				return fdtd::sources::chy * src.e[1][jpw_start_uni];
			});
		if (xpol)
			correct_plane(Hz, 1, jpw_near_h, subdomain.ipw_min, subdomain.ipw_max, subdomain.kpw_min, subdomain.kpw_max + 1, [&](int /*i*/, int /*k*/) {
				if (periodic())
					return -(fdtd::sources::chy * ref.e[0][jpw_start_uni]);
				return -(fdtd::sources::chy * src.e[0][jpw_start_uni]);
			});
		break;
	case Face::YFar:
		if (zpol)
			correct_plane(Hx, 1, subdomain.jpw_2, subdomain.ipw_min, subdomain.ipw_max + 1, subdomain.kpw_min, subdomain.kpw_max, [&](int, int) { return -(fdtd::sources::chy * src.e[1][jpw_end_uni]); });
		if (xpol)
			correct_plane(Hz, 1, subdomain.jpw_2, subdomain.ipw_min, subdomain.ipw_max, subdomain.kpw_min, subdomain.kpw_max + 1, [&](int, int) { return fdtd::sources::chy * src.e[0][jpw_end_uni]; });
		break;
	}
}

void tfsf_e(Face f)
{
	double c = 0;
	switch (f)
	{
	case Face::XNea: // Ey from Hz_inc, background of the y layer
	case Face::XFar:
	{
		const bool nea = f == Face::XNea;
		for (int j = subdomain.jpw_min; j < subdomain.jpw_max; ++j)
		{
			const int J = j + subdomain.j_for_uni;
			if (background_e_factor(fdtd::materials::background_y[J], c))
				correct_plane(Ey, 0, nea ? subdomain.ipw_1 : subdomain.ipw_2, j, j + 1, kpw_min_mod, subdomain.kpw_max + 1, [&](int, int) { return nea ? c * lat.inv_d[0] * src.h[0][J] : -(c * lat.inv_d[0] * src.h[0][J]); });
		}
		break;
	}
	case Face::ZNea: // Ey from Hx_inc, background of the y layer
	case Face::ZFar:
	{
		const bool nea = f == Face::ZNea;
		for (int j = subdomain.jpw_min; j < subdomain.jpw_max; ++j)
		{
			const int J = j + subdomain.j_for_uni;
			if (background_e_factor(fdtd::materials::background_y[J], c))
				correct_plane(Ey, 2, nea ? subdomain.kpw_1 : subdomain.kpw_2, ipw_min_mod, subdomain.ipw_max + 1, j, j + 1, [&](int, int) { return nea ? -(c * lat.inv_d[2] * src.h[1][J]) : c * lat.inv_d[2] * src.h[1][J]; });
		}
		break;
	}
	case Face::YNea: // first medium
		if (!background_e_factor(sp.medium[0], c))
			break;
		if (xpol)
			correct_plane(Ex, 1, subdomain.jpw_1, subdomain.ipw_min, subdomain.ipw_max, kpw_min_mod, subdomain.kpw_max + 1, [&](int /*i*/, int /*k*/) {
				if (periodic())
					return -(c * lat.inv_d[1] * ref.h[0][jpw_source_uni]);
				return -(c * lat.inv_d[1] * src.h[0][jpw_source_uni]);
			});
		if (zpol)
			correct_plane(Ez, 1, subdomain.jpw_1, ipw_min_mod, subdomain.ipw_max + 1, subdomain.kpw_min, subdomain.kpw_max, [&](int /*i*/, int /*k*/) {
				if (periodic())
					return c * lat.inv_d[1] * ref.h[1][jpw_source_uni];
				return c * lat.inv_d[1] * src.h[1][jpw_source_uni];
			});
		break;
	case Face::YFar: // second medium
		if (!background_e_factor(sp.medium[2], c))
			break;
		if (zpol)
			correct_plane(Ez, 1, subdomain.jpw_2, ipw_min_mod, subdomain.ipw_max + 1, subdomain.kpw_min, subdomain.kpw_max, [&](int, int) { return -(c * lat.inv_d[1] * src.h[1][jpw_end_uni]); });
		if (xpol)
			correct_plane(Ex, 1, subdomain.jpw_2, subdomain.ipw_min, subdomain.ipw_max, kpw_min_mod, subdomain.kpw_max + 1, [&](int, int) { return c * lat.inv_d[1] * src.h[0][jpw_end_uni]; });
		break;
	}
}

void tfsf_e_terms(int A, int i, int j, int k, double c, double &e_new)
{
	if (!sp.tfsf)
		return;
	if (A == 1)
	{
		// z planes, pol z, x planes, pol x: the incident H of the cell's y layer
		const bool on_j = j >= subdomain.jpw_min && j < subdomain.jpw_max;
		const int J = j + subdomain.j_for_uni;
		if (!lat.periodic_z && zpol && on_j && i >= ipw_min_mod && i <= subdomain.ipw_max)
		{
			if (subdomain.pw_nea_z && k == subdomain.kpw_1)
				e_new -= c * lat.inv_d[2] * src.h[1][J];
			if (subdomain.pw_far_z && k == subdomain.kpw_2)
				e_new += c * lat.inv_d[2] * src.h[1][J];
		}
		if (!lat.periodic_x && xpol && on_j && k >= kpw_min_mod && k <= subdomain.kpw_max)
		{
			if (subdomain.pw_nea_x && i == subdomain.ipw_1)
				e_new += c * lat.inv_d[0] * src.h[0][J];
			if (subdomain.pw_far_x && i == subdomain.ipw_2)
				e_new -= c * lat.inv_d[0] * src.h[0][J];
		}
		return;
	}
	// Ex, pol x, from Hz, and Ez, pol z, from Hx, on the y planes; with periodic x and z the near plane
	// takes the reference line
	const bool x = A == 0;
	if (!(x ? xpol : zpol))
		return;
	const bool in = x ? (i >= subdomain.ipw_min && i < subdomain.ipw_max && k >= kpw_min_mod && k <= subdomain.kpw_max) : (i >= ipw_min_mod && i <= subdomain.ipw_max && k >= subdomain.kpw_min && k < subdomain.kpw_max);
	if (!in)
		return;
	const double s = x ? -1 : 1; // Ex: - near, + far; Ez: + near, - far
	const cvector<double> &h = x ? src.h[0] : src.h[1], &h_ref = x ? ref.h[0] : ref.h[1];
	if (subdomain.pw_nea_y && j == subdomain.jpw_1)
	{
		const double v = c * lat.inv_d[1] * (lat.periodic_x && lat.periodic_z ? h_ref[jpw_source_uni] : h[jpw_source_uni]);
		if (s > 0)
			e_new += v;
		else
			e_new -= v;
	}
	if (subdomain.pw_far_y && j == subdomain.jpw_2)
	{
		const double v = c * lat.inv_d[1] * h[jpw_end_uni];
		if (s > 0)
			e_new -= v;
		else
			e_new += v;
	}
}
} // namespace fdtd::sources
