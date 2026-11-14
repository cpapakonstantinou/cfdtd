/**
 * @file materials/stretch.h
 * @brief Field updates over stretches: runs of cells along one k row with one CPML region and one material.
 *
 * @verbatim
      row (i, j), k ->
      | outer (CPML) |  inner: material runs from matind  ...........  | outer (CPML) |
      [k0          ia)[ia    )[      )[           )[                 ib)[          k1)
                      m = 4   m = 100   m = 4         m = 100
   @endverbatim
 *
 * A scan splits each row with for_each_stretch and hands every stretch to an update callback with the
 * material and whether it is inner, no CPML stretching in the curl. The update chooses a model:
 *   - Dielectric: E = c1 E + c2 curl, branch-free;
 *   - Drude2CP: the Drude + critical-points ADE over the whole stretch, for any component: the
 *     component's arrays come in an AdeRow, the curl and the CPML / TF-SF terms as callbacks;
 *   - another dispersive model: chosen once for the stretch, materials::with_model, its point
 *     update, inc/materials/models/, called for each cell;
 *   - otherwise the scan's per-cell update, nonlinear dielectrics, mode source.
 *
 * Every function does the same arithmetic in the same order as the per-cell code it replaces, so the
 * results are bit-identical, MARCH= builds.
 */
#ifndef FDTD_MATERIALS_STRETCH_H
#define FDTD_MATERIALS_STRETCH_H

#include "materials/codes.h"
#include "core/views.h"
#include "materials/material_field.h"
#include "materials/component.h"
#include "core/carray.h"
#include "materials/state.h"
#include <algorithm>
#include <vector>

namespace fdtd::materials::stretch
{
using fdtd::grid::subdomain;

/**
 * @brief Calls run(ka, kb, inner, material) for the stretches of the row k in [k0, k1)
 * @param ia, ib  the inner cells [ia, ib), whose material is row[k - org]; may be empty
 * @param outer   material of the other cells
 */
template <class Run>
inline void for_each_stretch(int k0, int k1, int ia, int ib, const int *row, int org, int outer, Run &&run)
{
	ia = std::clamp(ia, k0, k1);
	ib = std::clamp(ib, ia, k1);
	if (k0 < ia)
		run(k0, ia, false, outer);
	for (int k = ia; k < ib;)
	{
		const int m = row[k - org];
		int e = k + 1;
		while (e < ib && row[e - org] == m)
			++e;
		run(k, e, true, m);
		k = e;
	}
	if (ib < k1)
		run(ib, k1, false, outer);
}

/// The kernels below are always_inline: they run once per cell and must fold into the scan that
/// calls them, GCC stopped inlining Drude2CPPoint into the larger e_scan template: 155 M calls.

/// one stretch of a row: cells [k0, k1) of material m, inner = no CPML stretching in the curl
struct Stretch
{
	int k0, k1, m;
	bool inner;
};

/// bumped after every meshing, Driver::build_geometry; stretch lists built for an older mesh are rebuilt
inline int mesh_generation = 0;

/**
 * @brief The stretches of every row of one scan, split once per mesh instead of once per step
 *
 * A scan visits its rows in a fixed order. On the first step after a meshing it calls add_row() for
 * each row, the same arguments as for_each_stretch, and from then on for_row() replays the stored
 * stretches, finish() after the last row. Splitting costs about 0.3 ns per cell, a fifth of a dielectric update.
 */
class StretchList
{
	std::vector<Stretch> _stretches;
	std::vector<int> _row_end; ///< row r is _stretches[_row_end[r - 1], _row_end[r])
	int _generation = -1;

public:
	/// true when the list matches the current mesh; otherwise it is cleared for add_row()
	bool ready()
	{
		if (_generation == mesh_generation)
			return true;
		_stretches.clear(), _row_end.clear();
		return false;
	}
	/// splits and stores the next row
	void add_row(int k0, int k1, int ia, int ib, const int *row, int org, int outer)
	{
		for_each_stretch(k0, k1, ia, ib, row, org, outer, [&](int a, int b, bool inner, int m) { _stretches.push_back({a, b, m, inner}); });
		_row_end.push_back(int(_stretches.size()));
	}
	/// marks the list as built for the current mesh, after the last add_row
	void finish() { _generation = mesh_generation; }
	/// calls run(k0, k1, inner, m) for the stretches of row r
	template <class Run>
	void for_row(int r, Run &&run) const
	{
		for (int s = r ? _row_end[r - 1] : 0; s < _row_end[r]; ++s)
			run(_stretches[s].k0, _stretches[s].k1, _stretches[s].inner, _stretches[s].m);
	}
};

/// the arrays of one E component's row (i, j) for a dispersive material
struct AdeRow
{
	double *e, *e_prev; ///< E and E one step back
	double *p, *p_prev; ///< Drude polarisation and one step back
	double *pc[NrMaxCP], *pc_prev[NrMaxCP]; ///< critical-point polarisations and one step back
};

/// the AdeRow of row (i, j) of material mm for the arrays of one component
inline AdeRow make_ade_row(view3<double> e, view3<double> e_prev, material_field<double, nr_disp_mat_d> &p, material_field<double, nr_disp_mat_d> &p_prev,
						 ctetrad<double> *pc, ctetrad<double> *pc_prev, int mm, int i, int j)
{
	const auto &mp = fdtd::materials::model_parameters();
	AdeRow r{e[i][j], e_prev[i][j], p[mm][i][j], p_prev[mm][i][j], {}, {}};
	for (int n = 0; n < mp.nr_cp[mm]; ++n)
		r.pc[n] = pc[mm][n][i][j], r.pc_prev[n] = pc_prev[mm][n][i][j];
	return r;
}

/**
 * @brief CPML auxiliary field over [k0, k1) of one row
 * @verbatim
      psi[k - off] = c1[s(k)] psi[k - off] + c2[s(k)] dF(k)
      f[k]         = f[k] + factor(k) psi[k - off]
   @endverbatim
 * Every CPML face of every component, E and H, linear and nonlinear, is this update; only the arrays,
 * the grading index s, the derivative dF and the factor differ, H: +-coeff_h; E: +-ce_2 of the
 * background, or of the cell's material with mode_source. A factor -c gives f - c psi exactly.
 */
template <class S, class F, class D>
inline void psi_update(double *f, double *psi, int off, int k0, int k1, const cvector<double> &c1, const cvector<double> &c2, S &&s, F &&factor, D &&dF)
{
	for (int k = k0; k < k1; ++k)
	{
		double &p = psi[k - off];
		p = c1[s(k)] * p + c2[s(k)] * dF(k);
		f[k] = f[k] + factor(k) * p;
	}
}

} // namespace fdtd::materials::stretch

#endif // FDTD_STRETCH_UPDATE_H
