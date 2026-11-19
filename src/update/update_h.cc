/**
 * @file update/update_h.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-01-19
 * 
 * @brief Magnetic Field update functions
 * 
 * @details Houses the scnaning functions which update the magnetic field
 * components. Though there is a wrapper function for all components,
 * each component is updated by its own function.
 * 
 * @copyright Copyright 2022
 */

#include "boundary/cpml.h"
#include "materials/stretch.h"
#include "update/update.h"
#include <utility>
#include <algorithm>
#include <vector>
#include "grid/lattice.h"
#include "grid/stretching.h"
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

namespace
{
const fdtd::update::Spec &up = fdtd::update::spec();              ///< what the updates include
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace

namespace
{
fdtd::boundary::Cpml &the_cpml = fdtd::boundary::cpml(); ///< cached: the accessor is out of line
} // namespace

/**
 * @brief Main update function for the Magnetic Field
 * 
 * 
 * @details Updates each component of the magnetic field by first
 * managing special cases, TFSF, subdomain.has_plane_wave, then calling each
 * individual update scan function.
 */
void hx_update_scan();
void hy_update_scan();
void hz_update_scan();

void fdtd::update::update_h()
{

	hx_update_scan();
	hy_update_scan();
	hz_update_scan();

	return;
}

// Nested functions for the magnetic field components

namespace
{
/// CPML stretch factors on the half grid, H positions, along one axis: inner value inside [min, max)
template <class F>
double half_factor(const F &table, int n, int min, int max, int cells, int pad_nea, int pad_far, int num_cpml)
{
	return (n >= min && n < max) ? table[num_cpml] : n < min ? table[n + pad_nea] : table[cells - n - 1 + pad_far];
}

/// h += coeff_h curl(k) over [k0, k1) of one row
template <class Curl>
inline void add_curl(double *h, int k0, int k1, Curl &&curl)
{
	for (int k = k0; k < k1; ++k)
		h[k] += up.coeff_h * curl(k);
}

/// a CPML face row of an H component: the shared update with the factor sign * coeff_h
template <class S, class D>
inline void psi_update(double *h, double *psi, int off, int k0, int k1, const fdtd::array1<double> &c1, const fdtd::array1<double> &c2, S &&s, double sign, D &&dE)
{
	const double sc = sign * up.coeff_h;
	fdtd::materials::stretch::psi_update(h, psi, off, k0, k1, c1, c2, std::forward<S>(s), [=](int) { return sc; }, std::forward<D>(dE));
}

/// the z factors of every k of the half grid, they vary along a row only in the z CPML layers
std::vector<double> half_factors_z()
{
	std::vector<double> kz(std::max(subdomain.cells_z, 0));
	for (int k = 0; k < subdomain.cells_z; ++k)
		kz[k] = half_factor(fdtd::grid::stretching().inv_kappa[2][0][1], k, subdomain.k_cpml_min, subdomain.k_cpml_max, subdomain.cells_z,
						   subdomain.cpml_pad_nea_z, subdomain.cpml_pad_far_z, lat.cpml[2]);
	return kz;
}

/// the H arrays of one component: the field, the two E components of its curl, and its CPML auxiliary fields
struct HSet
{
	view3<double> h; ///< Hx: (Ey, Ez), Hy: (Ez, Ex), Hz: (Ex, Ey)
	view3<double> e1; ///< Hx: (Ey, Ez), Hy: (Ez, Ex), Hz: (Ex, Ey)
	view3<double> e2; ///< Hx: (Ey, Ez), Hy: (Ez, Ex), Hz: (Ex, Ey)
	view3<double> psi[4];              ///< near/far CPML of the first, then the second transverse axis
};

/**
 * @brief Hx += coeff_h (KZ dEy/dz - KY dEz/dy), then the y and z CPML faces
 * @verbatim
      i in [0, cells_x], j in [0, cells_y), k in [0, cells_z)
      psi: y near (-), y far (-), z near (+), z far (+)
   @endverbatim
 */
void update_hx(const HSet &f)
{
	const int nk = subdomain.cells_z;
	const std::vector<double> kz = half_factors_z();
	for (int i = 0; i <= subdomain.cells_x; ++i)
		for (int j = 0; j < subdomain.cells_y; ++j)
		{
			const double KY = half_factor(fdtd::grid::stretching().inv_kappa[1][0][1], j, subdomain.j_cpml_min, subdomain.j_cpml_max, subdomain.cells_y,
										subdomain.cpml_pad_nea_y, subdomain.cpml_pad_far_y, lat.cpml[1]);
			const double *ey = f.e1[i][j];
			const double *ez0 = f.e2[i][j];
			const double *ez1 = f.e2[i][j + 1];
			add_curl(f.h[i][j], 0, nk, [&](int k) { return kz[k] * ((ey[k + 1] - ey[k]) * lat.inv_d[2]) - KY * ((ez1[k] - ez0[k]) * lat.inv_d[1]); });
		}
	auto dEz_y = [&](int i, int j) { return [=](int k) { return (f.e2[i][j + 1][k] - f.e2[i][j][k]) * lat.inv_d[1]; }; };
	auto dEy_z = [&](int i, int j) { return [=](int k) { return (f.e1[i][j][k + 1] - f.e1[i][j][k]) * lat.inv_d[2]; }; };
	for (int i = 0; i <= subdomain.cells_x; ++i)
	{
		if (subdomain.nea_cpml_y())
			for (int j = 0; j < subdomain.j_nea_cpml_end; ++j)
			{
				const int s = j + subdomain.cpml_pad_nea_y;
				psi_update(f.h[i][j], f.psi[0][i][j], 0, 0, nk, the_cpml.profile[1][0].b[1], the_cpml.profile[1][0].c[1], [=](int) { return s; }, -1, dEz_y(i, j));
			}
		if (subdomain.far_cpml_y())
			for (int j = subdomain.j_far_cpml_start; j < subdomain.cells_y; ++j)
			{
				const int s = subdomain.cells_y - j - 1 + subdomain.cpml_pad_far_y;
				psi_update(f.h[i][j], f.psi[1][i][j - subdomain.j_far_cpml_start], 0, 0, nk, the_cpml.profile[1][1].b[1], the_cpml.profile[1][1].c[1], [=](int) { return s; }, -1, dEz_y(i, j));
			}
		for (int j = 0; j < subdomain.cells_y; ++j)
		{
			if (!lat.periodic_z && subdomain.nea_cpml_z())
				psi_update(f.h[i][j], f.psi[2][i][j], 0, 0, subdomain.k_nea_cpml_end, the_cpml.profile[2][0].b[1], the_cpml.profile[2][0].c[1],
					[](int k) { return k + subdomain.cpml_pad_nea_z; }, +1, dEy_z(i, j));
			if (!lat.periodic_z && subdomain.far_cpml_z())
				psi_update(f.h[i][j], f.psi[3][i][j], subdomain.k_far_cpml_start, subdomain.k_far_cpml_start, nk, the_cpml.profile[2][0].b[1], the_cpml.profile[2][0].c[1],
					[](int k) { return subdomain.cells_z - k - 1 + subdomain.cpml_pad_far_z; }, +1, dEy_z(i, j));
		}
	}
}

/**
 * @brief Hy += coeff_h (KX dEz/dx - KZ dEx/dz), then the x and z CPML faces
 * @verbatim
      i in [0, cells_x), j in [0, cells_y], k in [0, cells_z)
      psi: x near (+), x far (+), z near (-), z far (-)
   @endverbatim
 */
void update_hy(const HSet &f)
{
	const int nk = subdomain.cells_z;
	const std::vector<double> kz = half_factors_z();
	for (int i = 0; i < subdomain.cells_x; ++i)
	{
		const double KX = half_factor(fdtd::grid::stretching().inv_kappa[0][0][1], i, subdomain.i_cpml_min, subdomain.i_cpml_max, subdomain.cells_x,
									subdomain.cpml_pad_nea_x, subdomain.cpml_pad_far_x, lat.cpml[0]);
		for (int j = 0; j <= subdomain.cells_y; ++j)
		{
			const double *ez0 = f.e1[i][j];
			const double *ez1 = f.e1[i + 1][j];
			const double *ex = f.e2[i][j];
			add_curl(f.h[i][j], 0, nk, [&](int k) { return KX * ((ez1[k] - ez0[k]) * lat.inv_d[0]) - kz[k] * ((ex[k + 1] - ex[k]) * lat.inv_d[2]); });
		}
	}
	auto dEz_x = [&](int i, int j) { return [=](int k) { return (f.e1[i + 1][j][k] - f.e1[i][j][k]) * lat.inv_d[0]; }; };
	auto dEx_z = [&](int i, int j) { return [=](int k) { return (f.e2[i][j][k + 1] - f.e2[i][j][k]) * lat.inv_d[2]; }; };
	if (!lat.periodic_x && subdomain.nea_cpml_x())
		for (int i = 0; i < subdomain.i_nea_cpml_end; ++i)
		{
			const int s = i + subdomain.cpml_pad_nea_x;
			for (int j = 0; j <= subdomain.cells_y; ++j)
				psi_update(f.h[i][j], f.psi[0][i][j], 0, 0, nk, the_cpml.profile[0][0].b[1], the_cpml.profile[0][0].c[1], [=](int) { return s; }, +1, dEz_x(i, j));
		}
	if (!lat.periodic_x && subdomain.far_cpml_x())
		for (int i = subdomain.i_far_cpml_start; i < subdomain.cells_x; ++i)
		{
			const int s = subdomain.cells_x - i - 1 + subdomain.cpml_pad_far_x;
			for (int j = 0; j <= subdomain.cells_y; ++j)
				psi_update(f.h[i][j], f.psi[1][i - subdomain.i_far_cpml_start][j], 0, 0, nk, the_cpml.profile[0][0].b[1], the_cpml.profile[0][0].c[1], [=](int) { return s; }, +1, dEz_x(i, j));
		}
	for (int i = 0; i < subdomain.cells_x; ++i)
		for (int j = 0; j <= subdomain.cells_y; ++j)
		{
			if (!lat.periodic_z && subdomain.nea_cpml_z())
				psi_update(f.h[i][j], f.psi[2][i][j], 0, 0, subdomain.k_nea_cpml_end, the_cpml.profile[2][0].b[1], the_cpml.profile[2][0].c[1],
					[](int k) { return k + subdomain.cpml_pad_nea_z; }, -1, dEx_z(i, j));
			if (!lat.periodic_z && subdomain.far_cpml_z())
				psi_update(f.h[i][j], f.psi[3][i][j], subdomain.k_far_cpml_start, subdomain.k_far_cpml_start, nk, the_cpml.profile[2][0].b[1], the_cpml.profile[2][0].c[1],
					[](int k) { return subdomain.cells_z - k - 1 + subdomain.cpml_pad_far_z; }, -1, dEx_z(i, j));
		}
}

/**
 * @brief Hz += coeff_h (KY dEx/dy - KX dEy/dx), then the x and y CPML faces
 * @verbatim
      i in [0, cells_x), j in [0, cells_y), k in [0, cells_z]
      psi: x near (-), x far (-), y near (+), y far (+)
   @endverbatim
 */
void update_hz(const HSet &f)
{
	const int nk = subdomain.cells_z + 1;
	for (int i = 0; i < subdomain.cells_x; ++i)
	{
		const double KX = half_factor(fdtd::grid::stretching().inv_kappa[0][0][1], i, subdomain.i_cpml_min, subdomain.i_cpml_max, subdomain.cells_x,
									subdomain.cpml_pad_nea_x, subdomain.cpml_pad_far_x, lat.cpml[0]);
		for (int j = 0; j < subdomain.cells_y; ++j)
		{
			const double KY = half_factor(fdtd::grid::stretching().inv_kappa[1][0][1], j, subdomain.j_cpml_min, subdomain.j_cpml_max, subdomain.cells_y,
										subdomain.cpml_pad_nea_y, subdomain.cpml_pad_far_y, lat.cpml[1]);
			const double *ex0 = f.e1[i][j];
			const double *ex1 = f.e1[i][j + 1];
			const double *ey0 = f.e2[i][j];
			const double *ey1 = f.e2[i + 1][j];
			add_curl(f.h[i][j], 0, nk, [&](int k) { return KY * ((ex1[k] - ex0[k]) * lat.inv_d[1]) - KX * ((ey1[k] - ey0[k]) * lat.inv_d[0]); });
		}
	}
	auto dEy_x = [&](int i, int j) { return [=](int k) { return (f.e2[i + 1][j][k] - f.e2[i][j][k]) * lat.inv_d[0]; }; };
	auto dEx_y = [&](int i, int j) { return [=](int k) { return (f.e1[i][j + 1][k] - f.e1[i][j][k]) * lat.inv_d[1]; }; };
	if (!lat.periodic_x && subdomain.nea_cpml_x())
		for (int i = 0; i < subdomain.i_nea_cpml_end; ++i)
		{
			const int s = i + subdomain.cpml_pad_nea_x;
			for (int j = 0; j < subdomain.cells_y; ++j)
				psi_update(f.h[i][j], f.psi[0][i][j], 0, 0, nk, the_cpml.profile[0][0].b[1], the_cpml.profile[0][0].c[1], [=](int) { return s; }, -1, dEy_x(i, j));
		}
	if (!lat.periodic_x && subdomain.far_cpml_x())
		for (int i = subdomain.i_far_cpml_start; i < subdomain.cells_x; ++i)
		{
			const int s = subdomain.cells_x - i - 1 + subdomain.cpml_pad_far_x;
			for (int j = 0; j < subdomain.cells_y; ++j)
				psi_update(f.h[i][j], f.psi[1][i - subdomain.i_far_cpml_start][j], 0, 0, nk, the_cpml.profile[0][0].b[1], the_cpml.profile[0][0].c[1], [=](int) { return s; }, -1, dEy_x(i, j));
		}
	for (int i = 0; i < subdomain.cells_x; ++i)
	{
		if (subdomain.nea_cpml_y())
			for (int j = 0; j < subdomain.j_nea_cpml_end; ++j)
			{
				const int s = j + subdomain.cpml_pad_nea_y;
				psi_update(f.h[i][j], f.psi[2][i][j], 0, 0, nk, the_cpml.profile[1][0].b[1], the_cpml.profile[1][0].c[1], [=](int) { return s; }, +1, dEx_y(i, j));
			}
		if (subdomain.far_cpml_y())
			for (int j = subdomain.j_far_cpml_start; j < subdomain.cells_y; ++j)
			{
				const int s = subdomain.cells_y - j - 1 + subdomain.cpml_pad_far_y;
				psi_update(f.h[i][j], f.psi[3][i][j - subdomain.j_far_cpml_start], 0, 0, nk, the_cpml.profile[1][1].b[1], the_cpml.profile[1][1].c[1], [=](int) { return s; }, +1, dEx_y(i, j));
			}
	}
}
} // namespace

// The field and, when kept separately, the nonlinear field go through the same update with their own arrays.

void hx_update_scan()
{
	update_hx({Hx, Ey, Ez, {the_cpml.psi[0][1][0][1][0], the_cpml.psi[0][1][0][1][1], the_cpml.psi[0][1][0][2][0], the_cpml.psi[0][1][0][2][1]}});
	if (up.harmonic)
		update_hx({hx_nl, ey_nl, ez_nl, {the_cpml.psi[1][1][0][1][0], the_cpml.psi[1][1][0][1][1], the_cpml.psi[1][1][0][2][0], the_cpml.psi[1][1][0][2][1]}});
}

void hy_update_scan()
{
	update_hy({Hy, Ez, Ex, {the_cpml.psi[0][1][1][0][0], the_cpml.psi[0][1][1][0][1], the_cpml.psi[0][1][1][2][0], the_cpml.psi[0][1][1][2][1]}});
	if (up.harmonic)
		update_hy({hy_nl, ez_nl, ex_nl, {the_cpml.psi[1][1][1][0][0], the_cpml.psi[1][1][1][0][1], the_cpml.psi[1][1][1][2][0], the_cpml.psi[1][1][1][2][1]}});
}

void hz_update_scan()
{
	update_hz({Hz, Ex, Ey, {the_cpml.psi[0][1][2][0][0], the_cpml.psi[0][1][2][0][1], the_cpml.psi[0][1][2][1][0], the_cpml.psi[0][1][2][1][1]}});
	if (up.harmonic)
		update_hz({hz_nl, ex_nl, ey_nl, {the_cpml.psi[1][1][2][0][0], the_cpml.psi[1][1][2][0][1], the_cpml.psi[1][1][2][1][0], the_cpml.psi[1][1][2][1][1]}});
}
