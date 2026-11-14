/**
 * @file sources/aux_grid.cc
 * @brief The 1D auxiliary grid of the TF/SF plane-wave source: allocation and the E and H updates
 *
 * One update serves both polarisations and both lines, before: four copies in the E scan and two in
 * the H scan. The polarisations differ only in sign, pol 1: curl Ez = -dHx/dy, Hx -= , written
 * as add<P>, which adds for pol 0 and subtracts for pol 1, so the arithmetic is the same as before.
 */
#include "sources/aux_grid.h"
#include "sources/state.h"

#include "boundary/cpml.h"
#include "materials/models.h"
#include "grid/lattice.h"
#include "sources/spec.h"
#include "grid/stretching.h"
using fdtd::materials::DispMateOffsetD;
using fdtd::materials::Drude;
using fdtd::materials::DrudeCP;
using fdtd::materials::Lorentz;
using fdtd::materials::ModLorentz;
using fdtd::materials::NRMAXMATERIALS;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::sources; // transitional: module state, sources/state.h

namespace
{
fdtd::boundary::Cpml &the_cpml = fdtd::boundary::cpml(); ///< cached: the accessor is out of line
} // namespace

namespace fdtd::sources
{
AuxGrid &aux_grid()
{
	static AuxGrid g;
	return g;
}

void AuxLine::allocate(int n_cells, int n_cpml)
{
	auto zeroed = [](int n) {
		cvector<double> v = vector1<double>(n);
		for (int j = 0; j < n; ++j)
			v[j] = 0.;
		return v;
	};
	for (int q = 0; q < 2; ++q)
	{
		e[q] = zeroed(n_cells + 1);
		e_prev[q] = zeroed(n_cells + 1);
		h[q] = zeroed(n_cells + 1);
		p[q] = zeroed(n_cells + 1);
		p_prev[q] = zeroed(n_cells + 1);
		for (int n = 0; n < materials::NrMaxCP; ++n)
		{
			pc[q][n] = zeroed(n_cells + 1);
			pc_prev[q][n] = zeroed(n_cells + 1);
		}
		for (int side = 0; side < 2; ++side)
		{
			psi_e[q][side] = zeroed(n_cpml + 1);
			psi_h[q][side] = zeroed(n_cpml + 1);
		}
	}
}

void AuxLine::release()
{
	for (int q = 0; q < 2; ++q)
	{
		for (auto *v : {&e[q], &e_prev[q], &h[q], &p[q], &p_prev[q], &psi_e[q][0], &psi_e[q][1], &psi_h[q][0], &psi_h[q][1]})
			v->reset();
		for (int n = 0; n < materials::NrMaxCP; ++n)
		{
			pc[q][n].reset();
			pc_prev[q][n].reset();
		}
	}
}

void AuxGrid::release()
{
	mat.reset();
	line[0].release();
	line[1].release();
}

namespace
{
/// x += v for polarisation 0, x -= v for polarisation 1
template <int P>
inline void add(double &x, double v)
{
	if constexpr (P == 0)
		x += v;
	else
		x -= v;
}

/**
 * @brief E update of polarisation P at cell j of a line
 * @param m   material code of the cell
 * @param d   dH/dy of the polarisation's H, Hz for Ex, Hx for Ez
 * @param in  between the two CPMLs, no stretching; else the curl is scaled by K
 * @verbatim
      curl = +-d (K +-d in the CPML);  then the model of m, as the 3D models:
        Drude + CP, modified Lorentz:  CPML psi inside the E_new expression, scaled by c0 / ce3
        Drude, Lorentz, dielectric:    E updated, then CPML psi scaled by ce2
        other dispersive models:       no update, the 1D grid does not carry them
   @endverbatim
 */
template <int P>
void e_cell(AuxLine &L, int j, int m, double d, double K, bool in)
{
	const auto &mp = fdtd::materials::model_parameters();
	cvector<double> &e = L.e[P], &e_prev = L.e_prev[P], &p = L.p[P], &p_prev = L.p_prev[P];
	cvector<double>(&pc)[materials::NrMaxCP] = L.pc[P];
	cvector<double>(&pc_prev)[materials::NrMaxCP] = L.pc_prev[P];
	cvector<double> &psi_n = L.psi_e[P][0], &psi_f = L.psi_e[P][1];
	const double rotdiff = P == 0 ? (in ? d : d * K) : (in ? -d : -d * K);
	// CPML psi of the two ends, added to e_new scaled by c
	auto psi_terms = [&](double c, double &e_new) {
		if (j <= lat.cpml[1])
		{
			psi_n[j] = the_cpml.profile[1][0].b[0][j] * psi_n[j] + the_cpml.profile[1][0].c[0][j] * d;
			add<P>(e_new, c * psi_n[j]);
		}
		if (j >= nr_cell_par_uni)
		{
			const int s = sp.one_dim_cells - j, V = j - nr_cell_par_uni;
			psi_f[V] = the_cpml.profile[1][1].b[0][s] * psi_f[V] + the_cpml.profile[1][1].c[0][s] * d;
			add<P>(e_new, c * psi_f[V]);
		}
	};
	if (!(m >= DispMateOffsetD && m < NRMAXMATERIALS))
	{
		e[j] *= mp.ce_1[m];
		e[j] += mp.ce_2[m] * rotdiff;
		psi_terms(mp.ce_2[m], e[j]);
		return;
	}
	// the 3D models' point updates on the 1D line, a row indexed by j
	const int mm = m - DispMateOffsetD;
	const int nr_cp = mp.nr_cp[mm];
	fdtd::materials::stretch::AdeRow r{e, e_prev, p, p_prev, {}, {}};
	for (int n = 0; n < materials::NrMaxCP; ++n)
		r.pc[n] = pc[n], r.pc_prev[n] = pc_prev[n];
	if (mm < DrudeCP)
		fdtd::materials::models::Drude2CP::point(fdtd::materials::coeffs<fdtd::materials::models::Drude2CP>(mm), mm, j, r, rotdiff, psi_terms);
	else if (mm < Drude)
	{
		const auto &cf = fdtd::materials::coeffs<fdtd::materials::models::Drude>(mm);
		fdtd::materials::models::Drude::point(cf, j, r, rotdiff);
		psi_terms(cf.ce2, e[j]);
	}
	else if (mm < Lorentz)
	{
		const auto &cf = fdtd::materials::coeffs<fdtd::materials::models::Lorentz>(mm);
		fdtd::materials::models::Lorentz::point(cf, mm, j, r, rotdiff);
		psi_terms(cf.ce2, e[j]);
	}
	else if (mm < ModLorentz)
		fdtd::materials::models::ModifiedLorentz::point(fdtd::materials::coeffs<fdtd::materials::models::ModifiedLorentz>(mm), nr_cp, j, r, rotdiff, psi_terms);
}

/// H update of polarisation P at cell j; d = dE/dy of the polarisation's E
template <int P>
void h_cell(AuxLine &L, int j, double d, double K)
{
	cvector<double> &h = L.h[P], &psi_n = L.psi_h[P][0], &psi_f = L.psi_h[P][1];
	add<P>(h[j], sp.coeff_h * K * d);
	if (j < lat.cpml[1])
	{
		psi_n[j] = the_cpml.profile[1][0].b[1][j] * psi_n[j] + the_cpml.profile[1][0].c[1][j] * d;
		add<P>(h[j], sp.coeff_h * psi_n[j]);
	}
	if (j >= nr_cell_par_uni)
	{
		const int V = j - nr_cell_par_uni, s = nr_cell_uni - j;
		psi_f[V] = the_cpml.profile[1][1].b[1][s] * psi_f[V] + the_cpml.profile[1][1].c[1][s] * d;
		add<P>(h[j], sp.coeff_h * psi_f[V]);
	}
}

/// the lines that run: the source line, and the reference line with periodic x and z
int nr_lines() { return lat.periodic_x && lat.periodic_z ? 2 : 1; }
} // namespace

void AuxGrid::update_e()
{
	for (int l = 0; l < nr_lines(); ++l)
	{
		AuxLine &L = line[l];
		for (int j = 1; j < sp.one_dim_cells; ++j)
		{
			const bool in = j >= j_min_cpml_uni && j <= j_max_cpml_uni;
			const double K = fdtd::grid::stretching().inv_kappa[1][0][0][in ? lat.cpml[1] : j < j_min_cpml_uni ? j : sp.one_dim_cells - j];
			const int m = l == 0 ? mat[j] : sp.medium[0];
			e_cell<0>(L, j, m, (L.h[0][j] - L.h[0][j - 1]) * lat.inv_d[1], K, in);
			e_cell<1>(L, j, m, (L.h[1][j] - L.h[1][j - 1]) * lat.inv_d[1], K, in);
		}
	}
}

void AuxGrid::update_h()
{
	for (int l = 0; l < nr_lines(); ++l)
	{
		AuxLine &L = line[l];
		for (int j = 0; j < sp.one_dim_cells; ++j)
		{
			const bool in = j >= j_min_cpml_uni && j < j_max_cpml_uni;
			const double K = fdtd::grid::stretching().inv_kappa[1][0][1][in ? lat.cpml[1] : j < j_min_cpml_uni ? j : nr_cell_uni - j];
			h_cell<0>(L, j, (L.e[0][j + 1] - L.e[0][j]) * lat.inv_d[1], K);
			h_cell<1>(L, j, (L.e[1][j + 1] - L.e[1][j]) * lat.inv_d[1], K);
		}
	}
}
} // namespace fdtd::sources
