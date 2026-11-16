/**
 * @file materials/nonlinear/tensor.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @date 2022-06-16
 *
 * @brief The function is for adding the nonlinear perturbation to 
 * the, linear, polarization and electic field in the case of SHG and THG simulations.  

 * @details The exact update equations cannot be obtained for 
 * dispersive materials in the case of the full tensorial consideration of the nonlinear response. 
 * Therefore, in the calculation of the nonlinear perturbation, i.e., the second or third order nonlinear polarizations 
 * the current step electric field is assumed to be approximately equal to the previous step one: E_current = E_prev.

 * @copyright Copyright 2022
 */

//the approximation is P_i^{n+1} \approx P_i^{n}

#include "materials/nonlinear/nonlinear.h"
#include "core/kernels.h"
#include "materials/stretch.h"
#include "materials/models.h"
#include <array>
#include <vector>
#include "materials/spec.h"
using fdtd::constants::EP0;
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::fields::ex_prev;
using fdtd::fields::ex_temp;
using fdtd::fields::ey_prev;
using fdtd::fields::ey_temp;
using fdtd::fields::ez_prev;
using fdtd::fields::ez_temp;
using fdtd::grid::subdomain;
using fdtd::materials::DispMateOffsetD;
using fdtd::materials::Drude;
using fdtd::materials::DrudeCP;
using fdtd::materials::Lorentz;
using fdtd::materials::ModLorentz;
using fdtd::materials::NRMAXMATERIALS;
using fdtd::materials::matind_x;
using fdtd::materials::matind_y;
using fdtd::materials::matind_z;
using fdtd::materials::nl_code;
using fdtd::materials::nr_disp_mat_d;

namespace
{
} // namespace

namespace
{
/// the nonlinear cells of one scan, i, j, k, u, v, w, model index, listed once per mesh
struct NonlinearCells
{
	std::vector<std::array<int, 7>> cells;
	int generation = -1;
	/// rebuilds the list with fill(add) when the mesh changed; add(i, j, k, u, v, w, mm)
	template <class Fill>
	const std::vector<std::array<int, 7>> &get(Fill &&fill)
	{
		if (generation != fdtd::materials::stretch::mesh_generation)
		{
			cells.clear();
			fill([&](int i, int j, int k, int u, int v, int w, int mm) { cells.push_back({i, j, k, u, v, w, mm}); });
			generation = fdtd::materials::stretch::mesh_generation;
		}
		return cells;
	}
};

/**
 * @brief Nonlinear polarisation of one component, subtracted into the dispersive state of the
 *        nonlinear materials, code shifted by nl_code, of the inner region
 * @verbatim
      Drude-type, index < Lorentz:          P[mm]      -= P_NL
      modified Lorentz, Lorentz .. ModLor:  P_cp[mm][0] -= P_NL / dt
      P_NL = polarisation_term (inc/materials/nonlinear/nonlinear.h)
   @endverbatim
 * The loops cover the component's owned cells [i0, i1) x [j0, j1) x [k0, k1); only the inner
 * non-CPML, cells take part, as before. The nonlinear cells are listed once per mesh.
 */
void nonlinear_polarisation_scan(char c, int i0, int i1, int j0, int j1, int k0, int k1, view3<int> mat,
						   material_field<double, nr_disp_mat_d> &p, ctetrad<double> *pc, NonlinearCells &list)
{
	const auto &cells = list.get([&](auto &&add) {
		for (int i = i0; i < i1; ++i)
			for (int j = j0; j < j1; ++j)
				for (int k = k0; k < k1; ++k)
				{
					const int u = i - subdomain.inner_org_x;
					const int v = j - subdomain.inner_org_y;
					const int w = k - subdomain.inner_org_z;
					const int m = mat[u][v][w];
					if (m < nl_code + DispMateOffsetD || m >= nl_code + NRMAXMATERIALS)
						continue;
					const int mm = m - DispMateOffsetD - nl_code;
					if (fdtd::materials::nonlinear::nl[mm] && mm < ModLorentz)
						add(i, j, k, u, v, w, mm);
				}
	});
	for (const auto &[i, j, k, u, v, w, mm] : cells)
	{
		const int mx = matind_x[u][v][w];
		const int my = matind_y[u][v][w];
		const int mz = matind_z[u][v][w];
		const double p_nl = fdtd::materials::nonlinear::polarisation_term(mm, i, j, k, c - 'x', mx, my, mz, false);
		if (mm < Lorentz)
			p[mm][i][j][k] -= p_nl;
		else
			pc[mm][0][i][j][k] -= p_nl / fdtd::materials::model_spec().dt;
	}
}
} // namespace

/// Nonlinear polarisation of Ex, Ey and Ez, polarization approach, before the E update
void fdtd::materials::nonlinear::update_polarisation()
{
	auto &ds = fdtd::materials::dispersive_state();
	const auto &S = subdomain;
	static NonlinearCells lx;
	static NonlinearCells ly;
	static NonlinearCells lz;
	// each component's owned range intersected with its inner, non-CPML, cells
	nonlinear_polarisation_scan('x', std::max(0, S.i_cpml_min), std::min(S.cells_x, S.i_cpml_max), std::max(1, S.j_cpml_min), std::min(S.owned_cells_y, S.j_cpml_max) + 1,
						  std::max(1, S.k_cpml_min), std::min(S.owned_cells_z, S.k_cpml_max) + 1, matind_x, ds.c[0].p, ds.c[0].pc, lx);
	nonlinear_polarisation_scan('y', std::max(1, S.i_cpml_min), std::min(S.owned_cells_x, S.i_cpml_max) + 1, std::max(0, S.j_cpml_min), std::min(S.cells_y, S.j_cpml_max),
						  std::max(1, S.k_cpml_min), std::min(S.owned_cells_z, S.k_cpml_max) + 1, matind_y, ds.c[1].p, ds.c[1].pc, ly);
	nonlinear_polarisation_scan('z', std::max(1, S.i_cpml_min), std::min(S.owned_cells_x, S.i_cpml_max) + 1, std::max(1, S.j_cpml_min), std::min(S.owned_cells_y, S.j_cpml_max) + 1,
						  std::max(0, S.k_cpml_min), std::min(S.cells_z, S.k_cpml_max), matind_z, ds.c[2].p, ds.c[2].pc, lz);
}

namespace
{
/// f(i, j, k, u, v, w) over the inner, non-CPML, cells of [0, cells], (u, v, w) relative to the inner origin
template <class F>
void inner_cells(F &&f)
{
	const auto &S = subdomain;
	for (int i = std::max(0, S.i_cpml_min); i < std::min(S.cells_x + 1, S.i_cpml_max); ++i)
		for (int j = std::max(0, S.j_cpml_min); j < std::min(S.cells_y + 1, S.j_cpml_max); ++j)
			for (int k = std::max(0, S.k_cpml_min); k < std::min(S.cells_z + 1, S.k_cpml_max); ++k)
				f(i, j, k, i - S.inner_org_x, j - S.inner_org_y, k - S.inner_org_z);
}

/// model index of a nonlinear material that takes the E-approach update, code shifted by nl_code,
/// nonlinear, below the modified Lorentz range, or -1
int nonlinear_index(int m)
{
	if (m < nl_code + DispMateOffsetD || m >= nl_code + NRMAXMATERIALS)
		return -1;
	const int mm = m - DispMateOffsetD - nl_code;
	return (fdtd::materials::nonlinear::nl[mm] && mm < ModLorentz) ? mm : -1;
}
} // namespace

/**
 * @brief Nonlinear term added in the E update, E approach, nonlinear.update != polarization
 * @verbatim
      1) E_temp = E - coeff dP_NL   for every component, before any E changes, dP_NL reads Ex, Ey, Ez
         coeff: Drude2CP c1, Drude / Lorentz ce2 / dt, ModifiedLorentz ce3 / dt
      2) E = E_temp
   @endverbatim
 * Only the inner cells of nonlinear materials take part, listed once per mesh. dP_NL is fdtd::materials::nonlinear::polarisation_term(..., change = true).
 */
void fdtd::materials::nonlinear::update_field()
{
	struct Component
	{
		char c;
		view3<int> mat;
		view3<double> e;
		view3<double> temp;
	};
	const Component comps[3] = {{'x', matind_x, Ex, ex_temp}, {'y', matind_y, Ey, ey_temp}, {'z', matind_z, Ez, ez_temp}};
	static NonlinearCells lists[3];
	for (int a = 0; a < 3; ++a)
	{
		const Component &c = comps[a];
		for (const auto &[i, j, k, u, v, w, mm] : lists[a].get([&](auto &&add) {
				 inner_cells([&](int i, int j, int k, int u, int v, int w) {
					 const int mm = nonlinear_index(c.mat[u][v][w]);
					 if (mm >= 0)
						 add(i, j, k, u, v, w, mm);
				 });
			 }))
		{
			const double coeff = mm < DrudeCP ? double(fdtd::materials::coeffs<fdtd::materials::models::Drude2CP>(mm).c1) : mm < Lorentz ? double((mm < Drude ? fdtd::materials::coeffs<fdtd::materials::models::Drude>(mm).ce2 : fdtd::materials::coeffs<fdtd::materials::models::Lorentz>(mm).ce2) / fdtd::materials::model_spec().dt) : double(fdtd::materials::coeffs<fdtd::materials::models::ModifiedLorentz>(mm).ce3 / fdtd::materials::model_spec().dt);
			const int mx = matind_x[u][v][w];
			const int my = matind_y[u][v][w];
			const int mz = matind_z[u][v][w];
			const double dp = fdtd::materials::nonlinear::polarisation_term(mm, i, j, k, c.c - 'x', mx, my, mz, true);
			c.temp[i][j][k] = c.e[i][j][k] - coeff * dp;
		}
	}
	for (int a = 0; a < 3; ++a)
		for (const auto &[i, j, k, u, v, w, mm] : lists[a].cells)
			comps[a].e[i][j][k] = comps[a].temp[i][j][k];
}

double fdtd::materials::nonlinear::polarisation_term(int mm, int i, int j, int k, int a, int mx, int my, int mz, bool change)
{
	const double e[3] = {Ex[i][j][k], Ey[i][j][k], Ez[i][j][k]};
	const double ep[3] = {change ? double(ex_prev[i][j][k]) : double(0), change ? double(ey_prev[i][j][k]) : double(0), change ? double(ez_prev[i][j][k]) : double(0)};
	const int mat[3] = {mx, my, mz};
	const int b = a == 0 ? 1 : 0;
	const int c = a == 2 ? 1 : 2; // the other two axes, b < c
	const bool full = fdtd::materials::model_spec().staircase && mx == my && mx == mz;
	const bool with_b = mat[b] == mat[a];
	const bool with_c = mat[c] == mat[a];
	auto ee = [&](int p, int q) { return e[p] * e[q]; };
	auto ee_prev = [&](int p, int q) { return ep[p] * ep[q]; };
	// second-order product E_p E_q (p <= q), or its change over the step
	auto d = [&](int p, int q) { return change ? ee(p, q) - ee_prev(p, q) : ee(p, q); };
	// third-order factor of axis q: E_q^2, or the change of E_a E_q^2
	auto t = [&](int q) { return change ? e[a] * ee(q, q) - ep[a] * ee_prev(q, q) : ee(q, q); };
	auto voigt = [](int p, int q) { return p == q ? p : 6 - p - q; }; // xx yy zz | yz xz xy
	const double *d2 = chi_2[mm] + 6 * a;
	const double *c3 = chi_3[mm] + 7 * a;
	double temp = 0;
	if (fdtd::materials::model_spec().thg)
	{
		const double s_b = c3[1] + c3[2] + c3[3];
		const double s_c = c3[4] + c3[5] + c3[6];
		if (full)
			temp = change ? (c3[0] * t(a) + s_b * t(b) + s_c * t(c)) : e[a] * (c3[0] * t(a) + s_b * t(b) + s_c * t(c));
		else
		{
			temp = change ? c3[0] * t(a) : e[a] * c3[0] * t(a);
			if (with_b)
				temp += change ? s_b * t(b) : e[a] * s_b * t(b);
			if (with_c)
				temp += change ? s_c * t(c) : e[a] * s_c * t(c);
		}
	}
	if (fdtd::materials::model_spec().shg)
	{
		if (full)
			temp = (d2[0] * d(0, 0) + d2[1] * d(1, 1) + d2[2] * d(2, 2) + 2 * d2[3] * d(1, 2) + 2 * d2[4] * d(0, 2) + 2 * d2[5] * d(0, 1));
		else
		{
			temp = d2[a] * d(a, a);
			if (with_b)
				temp += (d2[b] * d(b, b) + 2 * d2[voigt(a, b)] * d(std::min(a, b), std::max(a, b)));
			if (with_c)
				temp += (d2[c] * d(c, c) + 2 * d2[voigt(a, c)] * d(std::min(a, c), std::max(a, c)));
			if (with_b && with_c)
				temp += 2 * d2[voigt(b, c)] * d(b, c);
		}
	}
	return EP0 * temp;
}
