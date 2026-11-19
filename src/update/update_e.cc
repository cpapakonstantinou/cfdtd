/**
 * @file update/update_e.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-01-19
 *
 * @brief Electric Field update functions
 *
 * @details Houses the scnaning functions which update the electric field
 * components. Though there is a wrapper function for all components,
 * each component is updated by its own function.
 *
 * @copyright Copyright 2022
 */

#include "io/text.h"
#include "update/update.h"
#include "materials/nonlinear/nonlinear.h"
#include "boundary/cpml.h"
#include "materials/models.h"
#include "materials/stretch.h"
#include <utility>
#include "grid/state.h"
#include "grid/lattice.h"
#include "grid/stretching.h"
using fdtd::constants::ROOTPROCESS;
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::fields::Hx;
using fdtd::fields::Hy;
using fdtd::fields::Hz;
using fdtd::fields::ex_nl;
using fdtd::fields::ex_prev;
using fdtd::fields::ex_prev_2;
using fdtd::fields::ex_temp;
using fdtd::fields::ey_nl;
using fdtd::fields::ey_prev;
using fdtd::fields::ey_prev_2;
using fdtd::fields::ey_temp;
using fdtd::fields::ez_nl;
using fdtd::fields::ez_prev;
using fdtd::fields::ez_prev_2;
using fdtd::fields::ez_temp;
using fdtd::fields::hx_nl;
using fdtd::fields::hy_nl;
using fdtd::fields::hz_nl;
using fdtd::grid::subdomain;
using fdtd::materials::max_slots;
using fdtd::materials::matind_x;
using fdtd::materials::matind_x_full;
using fdtd::materials::matind_y;
using fdtd::materials::matind_y_full;
using fdtd::materials::matind_z;
using fdtd::materials::matind_z_full;
using fdtd::materials::background_y;

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
 * @brief One E step: the 1D source grid, the nonlinear polarisation, P_NL approach, the three E
 *        scans, and the field history of the magnetic Lorentz model
 */
void ex_update_scan();
void ey_update_scan();
void ez_update_scan();

void fdtd::update::update_e()
{
	int i;
	int j;
	int k;

	if ((up.shg || up.thg) && up.pnl_polarisation)
	{
		static bool logged = false; // once per run
		if (fdtd::grid::world_rank == ROOTPROCESS && !logged && (logged = true))
		{
			fdtd::io::log("update", "NONLINEARITY IS ADDED TO THE POLARIZATION UPDATE (APPROACH: 1=POLARIZATION UPDATE).\n");
		}
		fdtd::materials::nonlinear::update_polarisation();
	}

	ex_update_scan();
	ey_update_scan();
	ez_update_scan();

	if (up.magnetic_lorentz)
	{
		for (i = 0; i <= subdomain.cells_x; i++)
		{
			for (j = 0; j <= subdomain.cells_y; j++)
			{
				for (k = 0; k <= subdomain.cells_z; k++)
				{
					ez_prev_2[i][j][k] = ez_prev[i][j][k];
					ez_prev[i][j][k] = ez_temp[i][j][k];
					ez_temp[i][j][k] = Ez[i][j][k];

					ey_prev_2[i][j][k] = ey_prev[i][j][k];
					ey_prev[i][j][k] = ey_temp[i][j][k];
					ey_temp[i][j][k] = Ey[i][j][k];

					ex_prev_2[i][j][k] = ex_prev[i][j][k];
					ex_prev[i][j][k] = ex_temp[i][j][k];
					ex_temp[i][j][k] = Ex[i][j][k];
				}
			}
		}
	}

	return;
}

// Nested functions for the electric field components

namespace
{
/**
 * @brief What differs between the Ex, Ey and Ez scans: arrays, model entry points, boundary terms
 *
 * The curl of E_A is p - q with p = dH_{A+2}/d(A+1), q = dH_{A+1}/d(A+2), cyclic x, y, z. The
 * per-component model functions and the CPML/TF-SF terms take the two derivatives in their own
 * order (dH1, dH2); the traits put (p, q) in that order.
 */
template <int A>
struct EAxis;

template <>
struct EAxis<0>
{
	static view3<double> e() { return Ex; }
	static view3<double> e_nl() { return ex_nl; }
	static view3<double> h1() { return Hy; } ///< H_{A+1}
	static view3<double> h2() { return Hz; } ///< H_{A+2}
	static view3<double> h1_nl() { return hy_nl; }
	static view3<double> h2_nl() { return hz_nl; }
	static view3<int> mat() { return matind_x; }
	static view3<int> mat_full() { return matind_x_full; }
	/// CPML psi / TF-SF terms: whole rows on the y CPML faces or the TF/SF y planes, cells on the z faces
	static bool row_terms(int, int j) { return (subdomain.nea_cpml_y() && j <= subdomain.j_nea_cpml_end) || (subdomain.far_cpml_y() && j >= subdomain.j_far_cpml_start) || (up.tfsf && up.tfsf_xpol && (j == subdomain.jpw_1 || j == subdomain.jpw_2)); }
	static bool k_terms(int k) { return !lat.periodic_z && ((subdomain.nea_cpml_z() && k <= subdomain.k_nea_cpml_end) || (subdomain.far_cpml_z() && k >= subdomain.k_far_cpml_start)); }
	static void shg_dl(double r, int m, int i, int j, int k) { fdtd::materials::nonlinear::update_exact_shg(0, r, m, i, j, k); }
	static void thg_dl(double r, int m, int i, int j, int k) { fdtd::materials::nonlinear::update_exact_thg(0, r, m, i, j, k); }
	static void psi() { fdtd::boundary::update_psi_ex(); }
};

template <>
struct EAxis<1>
{
	static view3<double> e() { return Ey; }
	static view3<double> e_nl() { return ey_nl; }
	static view3<double> h1() { return Hz; }
	static view3<double> h2() { return Hx; }
	static view3<double> h1_nl() { return hz_nl; }
	static view3<double> h2_nl() { return hx_nl; }
	static view3<int> mat() { return matind_y; }
	static view3<int> mat_full() { return matind_y_full; }
	/// CPML psi / TF-SF terms: whole rows on the x CPML faces or the TF/SF x planes, cells on the z faces or the TF/SF z planes
	static bool row_terms(int i, int) { return (!lat.periodic_x && subdomain.nea_cpml_x() && i <= subdomain.i_nea_cpml_end) || (!lat.periodic_x && subdomain.far_cpml_x() && i >= subdomain.i_far_cpml_start) || (up.tfsf && up.tfsf_xpol && (i == subdomain.ipw_1 || i == subdomain.ipw_2)); }
	static bool k_terms(int k) { return (!lat.periodic_z && ((subdomain.nea_cpml_z() && k <= subdomain.k_nea_cpml_end) || (subdomain.far_cpml_z() && k >= subdomain.k_far_cpml_start))) || (up.tfsf && up.tfsf_zpol && (k == subdomain.kpw_1 || k == subdomain.kpw_2)); }
	static void shg_dl(double r, int m, int i, int j, int k) { fdtd::materials::nonlinear::update_exact_shg(1, r, m, i, j, k); }
	static void thg_dl(double r, int m, int i, int j, int k) { fdtd::materials::nonlinear::update_exact_thg(1, r, m, i, j, k); }
	static void psi() { fdtd::boundary::update_psi_ey(); }
};

template <>
struct EAxis<2>
{
	static view3<double> e() { return Ez; }
	static view3<double> e_nl() { return ez_nl; }
	static view3<double> h1() { return Hx; }
	static view3<double> h2() { return Hy; }
	static view3<double> h1_nl() { return hx_nl; }
	static view3<double> h2_nl() { return hy_nl; }
	static view3<int> mat() { return matind_z; }
	static view3<int> mat_full() { return matind_z_full; }
	/// CPML psi / TF-SF terms: whole rows on the x and y CPML faces or the TF/SF y planes, no k-dependent terms
	static bool row_terms(int i, int j)
	{
		return (subdomain.nea_cpml_y() && j <= subdomain.j_nea_cpml_end) || (subdomain.far_cpml_y() && j >= subdomain.j_far_cpml_start) ||
			   (!lat.periodic_x && subdomain.nea_cpml_x() && i <= subdomain.i_nea_cpml_end) || (!lat.periodic_x && subdomain.far_cpml_x() && i >= subdomain.i_far_cpml_start) ||
			   (up.tfsf && up.tfsf_zpol && (j == subdomain.jpw_1 || j == subdomain.jpw_2));
	}
	static bool k_terms(int) { return false; }
	static void shg_dl(double r, int m, int i, int j, int k) { fdtd::materials::nonlinear::update_exact_shg(2, r, m, i, j, k); }
	static void thg_dl(double r, int m, int i, int j, int k) { fdtd::materials::nonlinear::update_exact_thg(2, r, m, i, j, k); }
	static void psi() { fdtd::boundary::update_psi_ez(); }
};

/// 1 / kappa of the CPML at index n along axis b, the interface value inside the domain
double cpml_factor(int b, int n, bool in)
{
	const auto &S = subdomain;
	const int cells[3] = {S.cells_x, S.cells_y, S.cells_z};
	const int min[3] = {S.i_cpml_min, S.j_cpml_min, S.k_cpml_min};
	const int pad_nea[3] = {S.cpml_pad_nea_x, S.cpml_pad_nea_y, S.cpml_pad_nea_z};
	const int pad_far[3] = {S.cpml_pad_far_x, S.cpml_pad_far_y, S.cpml_pad_far_z};
	const int cpml[3] = {lat.cpml[0], lat.cpml[1], lat.cpml[2]};
	const fdtd::array1<double> &t = b == 0 ? fdtd::grid::stretching().inv_kappa[0][0][0] : b == 1 ? fdtd::grid::stretching().inv_kappa[1][0][0] : fdtd::grid::stretching().inv_kappa[2][0][0];
	return t[in ? cpml[b] : n < min[b] ? n + pad_nea[b] : cells[b] - n + pad_far[b]];
}
/// backward difference of h along axis b at (i, j, k), times 1 / d_b
template <int b>
inline double diff(const view3<double> &h, int i, int j, int k)
{
	if constexpr (b == 0)
		return (h[i][j][k] - h[i - 1][j][k]) * lat.inv_d[0];
	else if constexpr (b == 1)
		return (h[i][j][k] - h[i][j - 1][k]) * lat.inv_d[1];
	else
		return (h[i][j][k] - h[i][j][k - 1]) * lat.inv_d[2];
}
} // namespace

/**
 * @brief Update of the E component A (0 x, 1 y, 2 z) on every owned cell, row by row in stretches
 * @verbatim
      rows (i, j), cells k; along A the loop covers 0 .. cells and the inner test is exclusive
      n < max, across A it covers 1 .. owned and the inner test is inclusive, n <= max
      curl  = p - q, inner cells,  K_{A+1} p - K_{A+2} q, CPML,
              p = dH_{A+2}/d(A+1),  q = dH_{A+1}/d(A+2),  K = 1/kappa of the CPML
      stretches, inc/update/stretch.h: a model with a row loop, dielectric, local Drude + 2CP, modified
      Lorentz, Drude, Lorentz, takes the whole stretch, every other material the per-cell update;
      then the CPML psi of the E faces
   @endverbatim
 * The inner cells take their material from the component's material map, the CPML cells the
 * background of their y layer.
 */
template <int A>
void e_scan()
{
	const auto &mp = fdtd::materials::model_parameters();
	using T = EAxis<A>;
	static fdtd::materials::stretch::StretchList list; // rows in the order of the loops below
	const bool replay = list.ready();
	const auto &S = subdomain;
	const view3<double> E = T::e();
	const view3<double> E_nl = T::e_nl();
	const view3<double> h1 = T::h1();
	const view3<double> h2 = T::h2();
	const view3<double> h1_nl = T::h1_nl();
	const view3<double> h2_nl = T::h2_nl();
	constexpr int B1 = (A + 1) % 3;
	constexpr int B2 = (A + 2) % 3; // the two derivative axes
	int row = 0;
	const int i_first = A == 0 ? 0 : 1;
	const int i_last = A == 0 ? S.cells_x : S.owned_cells_x;
	const int j_first = A == 1 ? 0 : 1;
	const int j_last = A == 1 ? S.cells_y : S.owned_cells_y;
	for (int i = i_first; i <= i_last; ++i)
	{
		const bool u_in = A == 0 ? (i >= S.i_cpml_min && i < S.i_cpml_max) : (i >= S.i_cpml_min && i <= S.i_cpml_max);
		const int u = i - S.inner_org_x;
		for (int j = j_first; j <= j_last; ++j)
		{
			const int J2 = j + subdomain.output_start[1];
			const bool v_in = A == 1 ? (j >= S.j_cpml_min && j < S.j_cpml_max) : (j >= S.j_cpml_min && j <= S.j_cpml_max);
			const int m_outer = v_in ? background_y[j + subdomain.j_for_uni] : j < S.j_cpml_min ? up.medium_near_y : up.medium_far_y;
			// 1/kappa along the derivative axes: x and y once per row, z per cell
			const double Kx = A != 0 ? cpml_factor(0, i, u_in) : 0;
			const double Ky = A != 1 ? cpml_factor(1, j, v_in) : 0;
			auto Kz = [&](int k) {
				return fdtd::grid::stretching().inv_kappa[2][0][0][(k >= S.k_cpml_min && k <= S.k_cpml_max) ? lat.cpml[2] : k < S.k_cpml_min ? k + S.cpml_pad_nea_z : S.cells_z - k + S.cpml_pad_far_z];
			};
			auto K1 = [&](int k) { if constexpr (B1 == 0) return double(Kx); else if constexpr (B1 == 1) return double(Ky); else return double(Kz(k)); };
			auto K2 = [&](int k) { if constexpr (B2 == 0) return double(Kx); else if constexpr (B2 == 1) return double(Ky); else return double(Kz(k)); };
			// (p, q) of the field or of the harmonic field
			auto derivatives = [&](int k, bool nl) {
				return nl ? std::pair<double, double>{diff<B1>(h2_nl, i, j, k), diff<B2>(h1_nl, i, j, k)}
						  : std::pair<double, double>{diff<B1>(h2, i, j, k), diff<B2>(h1, i, j, k)};
			};
			auto combine = [&](int k, bool inner, double p, double q) { return inner ? p - q : K1(k) * p - K2(k) * q; };
			// one cell of a dispersive model; the CPML psi and TF/SF terms are the boundary callback
			auto model_cell = [&]<class M>(int mm, int k, double rotdiff, double rotdiff_NL, double p, double q) {
				M::template update<A>(fdtd::materials::models::Cell{mm, i, j, k, J2, rotdiff, rotdiff_NL},
									  [&](double c, double &e_new) { fdtd::boundary::cpml_tfsf_terms<A>(i, j, k, c, p, q, e_new); });
			};
			// the per-cell update, all material models, CPML psi and TF/SF terms
			auto cell = [&](int k, bool inner, int m) {
				const auto [p, q] = derivatives(k, false);
				const double rotdiff = combine(k, inner, p, q);
				double rotdiff_NL = 0;
				if (up.harmonic)
				{
					const auto [n1, n2] = derivatives(k, true);
					rotdiff_NL = combine(k, inner, n1, n2);
				}
				if (fdtd::materials::dispersive(m))
					fdtd::materials::with_model(m, [&]<class M>(M) { model_cell.template operator()<M>(m, k, rotdiff, rotdiff_NL, p, q); });
				else
				{
					if (fdtd::materials::nonlinear::nl[m] && (fdtd::materials::nonlinear::iso_chi_2[m] || fdtd::materials::nonlinear::iso_chi_3[m]))
					{
						if (up.shg)
							T::shg_dl(rotdiff, m, i, j, k);
						if (up.thg)
							T::thg_dl(rotdiff, m, i, j, k);
					}
					else if (!up.cellwise_materials)
					{
						E[i][j][k] *= mp.ce_1[m];
						E[i][j][k] += mp.ce_2[m] * rotdiff;
						if (up.harmonic)
						{
							E_nl[i][j][k] *= mp.ce_1[m];
							E_nl[i][j][k] += mp.ce_2[m] * rotdiff_NL;
						}
					}
					else
					{
						const int mc = T::mat_full()[i][j][k];
						E[i][j][k] *= mp.ce_1[mc];
						E[i][j][k] += mp.ce_2[mc] * rotdiff;
					}
				}
			};
			// one stretch: a branch-free model when one applies, else the per-cell update
			auto update = [&](int k0, int k1, bool inner, int m) {
				auto curl = [&](int k) { return combine(k, inner, diff<B1>(h2, i, j, k), diff<B2>(h1, i, j, k)); }; // the hot loop: no nl branch
				using fdtd::materials::models::Dielectric;
				using fdtd::materials::models::Drude2CP;
				// the CPML psi and TF/SF terms of cell k: only on the rows and cells that have them
				const bool row_terms = T::row_terms(i, j);
				auto boundary = [&](int k, double c, double &e_new) {
					if (row_terms || T::k_terms(k))
						fdtd::boundary::cpml_tfsf_terms<A>(i, j, k, c, diff<B1>(h2, i, j, k), diff<B2>(h1, i, j, k), e_new);
				};
				if (Dielectric::applies(m))
					Dielectric::stretch(E[i][j], k0, k1, mp.ce_1[m], mp.ce_2[m], curl);
				else if (const int mm = Drude2CP::stretch_index(m); mm >= 0)
					Drude2CP::stretch<A>(mm, i, j, k0, k1, curl, boundary);
				else if (fdtd::materials::dispersive(m))
				{
					// another dispersive model, chosen once for the stretch: its row loop when it has one
					// with the boundary callback: modified Lorentz; without: Drude, Lorentz, else the
					// per-cell update. Drude + 2CP arrives here only when its row path is excluded
					// nonlocal, electron-density DFT; Drude2CP::stretch_index.
					const int mm = m;
					fdtd::materials::with_model(mm, [&]<class M>(M) {
						if constexpr (!std::is_same_v<M, Drude2CP> && requires { M::template stretch<A>(mm, i, j, k0, k1, curl, boundary); })
						{
							M::template stretch<A>(mm, i, j, k0, k1, curl, boundary);
							return;
						}
						else if constexpr (requires { M::template stretch<A>(mm, i, j, k0, k1, curl); })
						{
							M::template stretch<A>(mm, i, j, k0, k1, curl);
							return;
						}
						for (int k = k0; k < k1; ++k)
						{
							const auto [p, q] = derivatives(k, false);
							const double rotdiff = combine(k, inner, p, q);
							double rotdiff_NL = 0;
							if (up.harmonic)
							{
								const auto [n1, n2] = derivatives(k, true);
								rotdiff_NL = combine(k, inner, n1, n2);
							}
							model_cell.template operator()<M>(mm, k, rotdiff, rotdiff_NL, p, q);
						}
					});
				}
				else
					for (int k = k0; k < k1; ++k)
						cell(k, inner, m);
			};
			// inner cells along k: [k_cpml_min, k_cpml_max] across z, [k_cpml_min, k_cpml_max) along z (Ez)
			const int k_first = A == 2 ? 0 : 1;
			const int k_end = A == 2 ? S.cells_z + 1 : S.owned_cells_z + 1;
			const int k0 = S.k_cpml_min;
			const int k1 = A == 2 ? S.k_cpml_max : S.k_cpml_max + 1;
			const bool in = u_in && v_in && S.inner_dim_z > 0; // a rank that is all CPML along z has no matind
			if (!replay)
				list.add_row(k_first, k_end, in ? k0 : k_first, in ? k1 : k_first, in ? T::mat()[u][j - S.inner_org_y] : nullptr, S.inner_org_z, m_outer);
			list.for_row(row++, update);
		}
	}
	if (!replay)
		list.finish();
	T::psi();
}

void ex_update_scan() { e_scan<0>(); }
void ey_update_scan() { e_scan<1>(); }
void ez_update_scan() { e_scan<2>(); }
