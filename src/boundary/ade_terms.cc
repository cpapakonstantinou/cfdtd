/**
 * @file boundary/ade_terms.cc
 * @brief The CPML and TF/SF terms of the E update of a dispersive cell, one template for Ex, Ey, Ez.
 *
 * Every ADE model that takes the boundary callback, Drude + critical points, modified Lorentz,
 * nonlinear Drude-Lorentz, hydrodynamic, adds these terms to its new E value, scaled by its own
 * coefficient. The CPML psi of the two derivative axes is one function, psi_update<A, D, sign>; the
 * TF/SF terms are the incident field the application hands over, boundary/incident_field.h.
 */
#include "boundary/cpml.h"
#include "boundary/incident_field.h"
#include "grid/lattice.h"
using fdtd::boundary::index::far;
using fdtd::grid::subdomain;

namespace
{
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace

namespace fdtd::boundary
{
namespace
{
Cpml &the_cpml = cpml(); ///< this rank's CPML, cpml() is out of line; this is read per cell
IncidentETerms incident_e = nullptr; ///< set_incident_e_terms

/// x += v (S > 0) or x -= v (S < 0)
template <int S>
inline void add(double &x, double v)
{
	if constexpr (S > 0)
		x += v;
	else
		x -= v;
}

/**
 * @brief CPML psi of component A along derivative axis D at (i, j, k), both faces, added with sign S
 * @verbatim
      near face, n <= nea_end:          s = n + pad_nea          psi[ijk]
      far face,  far_start <= n, < cells: s = cells - n + pad_far  psi[ijk shifted by -far_start along D]
      psi = b[s] psi + c[s] d;   e_new += S c psi
   @endverbatim
 * The x and z faces take the near-face profile on both faces, the y faces their own, as before.
 */
template <int A, int D, int S>
inline void psi_update(int i, int j, int k, double coeff, double d, double &e_new)
{
	const auto &sd = subdomain;
	if (D == 0 && lat.periodic_x)
		return;
	if (D == 2 && lat.periodic_z)
		return;
	const int n = D == 0 ? i : D == 1 ? j : k;
	const bool nea = D == 0 ? subdomain.nea_cpml_x() : D == 1 ? subdomain.nea_cpml_y() : subdomain.nea_cpml_z();
	const bool far = D == 0 ? subdomain.far_cpml_x() : D == 1 ? subdomain.far_cpml_y() : subdomain.far_cpml_z();
	const int nea_end = D == 0 ? sd.i_nea_cpml_end : D == 1 ? sd.j_nea_cpml_end : sd.k_nea_cpml_end;
	const int far_start = D == 0 ? sd.i_far_cpml_start : D == 1 ? sd.j_far_cpml_start : sd.k_far_cpml_start;
	const int pad_nea = D == 0 ? sd.cpml_pad_nea_x : D == 1 ? sd.cpml_pad_nea_y : sd.cpml_pad_nea_z;
	const int pad_far = D == 0 ? sd.cpml_pad_far_x : D == 1 ? sd.cpml_pad_far_y : sd.cpml_pad_far_z;
	const int cells = D == 0 ? sd.cells_x : D == 1 ? sd.cells_y : sd.cells_z;
	Cpml &cp = the_cpml;
	if (nea && n <= nea_end)
	{
		const int s = n + pad_nea;
		const Profile &P = cp.profile[D][0];
		double &psi = cp.psi[index::linear][index::e][A][D][0][i][j][k];
		psi = P.b[index::integer][s] * psi + P.c[index::integer][s] * d;
		add<S>(e_new, coeff * psi);
	}
	if (far && n >= far_start && ((pad_far > 0 && n <= cells) || (pad_far == 0 && n < cells)) && n != 0)
	{
		const int s = cells - n + pad_far;
		const Profile &P = cp.profile[D][D == 1 ? 1 : 0];
		int x[3] = {i, j, k};
		x[D] -= far_start;
		double &psi = cp.psi[index::linear][index::e][A][D][1][x[0]][x[1]][x[2]];
		psi = P.b[index::integer][s] * psi + P.c[index::integer][s] * d;
		add<S>(e_new, coeff * psi);
	}
}
} // namespace

template <int A>
void cpml_tfsf_terms(int i, int j, int k, double c, double p, double q, double &e_new)
{
	constexpr int B1 = (A + 1) % 3;
	constexpr int B2 = (A + 2) % 3;
	// curl = p - q: the B1 faces add, the B2 faces subtract; summed y first (Ex: y, z), else x first
	// (Ey: x, z; Ez: x is B1, y is B2 -> y, x), the order of the previous code
	if constexpr (A == 0)
	{
		psi_update<A, B1, +1>(i, j, k, c, p, e_new);
		psi_update<A, B2, -1>(i, j, k, c, q, e_new);
	}
	else
	{
		psi_update<A, B2, -1>(i, j, k, c, q, e_new);
		psi_update<A, B1, +1>(i, j, k, c, p, e_new);
	}
	if (incident_e)
		incident_e(A, i, j, k, c, e_new);
}

void set_incident_e_terms(IncidentETerms terms) { incident_e = terms; }

template void cpml_tfsf_terms<0>(int, int, int, double, double, double, double &);
template void cpml_tfsf_terms<1>(int, int, int, double, double, double, double &);
template void cpml_tfsf_terms<2>(int, int, int, double, double, double, double &);
} // namespace fdtd::boundary
