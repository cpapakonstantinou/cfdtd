/**
 * @file materials/models/drude_2cp.h
 * @brief Drude + critical points model, local and nonlocal, model drude_cp
 *
 * point() is the local arithmetic, which the scans run over whole stretches, see stretch(), and the 1D
 * auxiliary grid per cell; update() adds the nonlocal / GNOR variant and the electron-density sum.
 */
#ifndef FDTD_MATERIALS_MODELS_DRUDE_2CP_H
#define FDTD_MATERIALS_MODELS_DRUDE_2CP_H

#include "materials/component.h"
#include "materials/state.h"
#include "materials/models/cell.h"
#include "materials/coefficients.h"
#include <cmath>
#include "materials/stretch.h"
#include "materials/spec.h"
#include "grid/lattice.h"

namespace fdtd::materials::models
{
/**
 * @brief Drude + critical points ADE
 * @verbatim
      E_new  = c0 curl + c1 (-a0 P_prev + (1 - a1) P) + c2 E_prev + c3 E          [local]
               nonlocal: P from p_nonlocal, minus c0 non_loc grad_div / DT, GNOR correction
             + sum_n c1 (-b0_n Pn_prev + (1 - b1_n) Pn) + CPML/TF-SF terms * c0
      P_new  = a0 P_prev + a1 P + a2 (E_prev + 2 E + E_new) (+ non_loc grad_div, GNOR)
      Pn_new = b0_n Pn_prev + b1_n Pn + b2_n E_prev + b3_n E + b4_n E_new
   @endverbatim
 */
struct Drude2CP
{
	/// update coefficients of one material: the poles and the Drude pole (a0, a1, a2)
	struct Coeffs : Poles
	{
		double a0; ///< P_new = a0 P_prev + a1 P + a2 (E_prev + 2 E + E_new)
		double a1; ///< P_new = a0 P_prev + a1 P + a2 (E_prev + 2 E + E_new)
		double a2; ///< P_new = a0 P_prev + a1 P + a2 (E_prev + 2 E + E_new)
	};

	/// coefficients of material slot mm for time step dt;
	/// c0..c3 take the sums of b2, b3, b4 over this material's poles
	static Coeffs compute(int mm, double dt)
	{
		const auto &mp = fdtd::materials::model_parameters();
		double acc[3] = {0, 0, 0};
		using namespace fdtd::constants; // EP0, MU0, ...
		const int n = mm;
		double cp1;
		double cp2;
		double cp3;
		double cp4;
		double cp5;
		double cp6;
		double cp7;
		double cp8;
		double cp9;
		double cp10;
		Coeffs c{};
		const double alpha = mp.ome_plasma[n] * mp.ome_plasma[n] * dt * dt;
		const double beta = mp.gam_plasma[n] * dt;
		c.a0 = (beta - 2) / (beta + 2);
		c.a1 = 4 / (beta + 2);
		c.a2 = EP0 * alpha / (2 * (beta + 2));
		for (int s = 0; s < mp.nr_cp[n]; ++s)
		{
			cp1 = mp.big_omega_cp[n][s] * mp.big_omega_cp[n][s] * dt * dt;
			cp2 = mp.big_gamma_cp[n][s] * dt;
			cp3 = EP0 * mp.a_coeff_cp[n][s] * mp.big_omega_cp[n][s] * dt;
			cp4 = (mp.big_omega_cp[n][s] * mp.big_omega_cp[n][s] + mp.big_gamma_cp[n][s] * mp.big_gamma_cp[n][s]) * dt * dt;
			cp5 = mp.big_omega_cp[n][s] * dt;
			cp6 = cp1 + (2 + cp2) * (2 + cp2);
			cp7 = cp5 * std::cos(mp.phi_cp[n][s]);
			cp8 = std::sin(mp.phi_cp[n][s]);

			c.b0[s] = -(cp1 + (2 - cp2) * (2 - cp2)) / cp6;
			c.b1[s] = 2 * (4 - cp4) / cp6;
			c.b2[s] = 2 * cp3 * (cp7 + (2 - cp2) * cp8) / cp6;
			c.b3[s] = 4 * cp3 * (cp7 - cp2 * cp8) / cp6;
			c.b4[s] = 2 * cp3 * (cp7 - (2 + cp2) * cp8) / cp6;
			acc[0] += c.b2[s];
			acc[1] += c.b3[s];
			acc[2] += c.b4[s];
		}
		cp9 = mp.sig_m[n] * dt / 2 + 2 * c.a2 + acc[1] - EP0 * mp.epre_inf[n];
		cp10 = mp.sig_m[n] * dt / 2 + c.a2 + acc[2] + EP0 * mp.epre_inf[n];
		c.c0 = dt / cp10;
		c.c1 = 1 / cp10;
		c.c2 = -(c.a2 + acc[0]) / cp10;
		c.c3 = -cp9 / cp10;
		return c;
	}

	/**
	 * @brief grad(div P) of component A at the padded index x, from the nonlocal polarisations
	 * @verbatim
	      d_A d_A P_A + sum_{B != A} d_A d_B P_B, backward-forward differences on the staggered grid
	   @endverbatim
	 * The two cross terms are summed y first, or x first for A = y: the order of the previous build.
	 */
	template <int A>
	static double grad_div(const int (&x)[3])
	{
		constexpr int F = A == 1 ? 0 : 1;
		constexpr int S = 3 - A - F;
		const view3<double> pa = Component<A>::p_nonlocal_2();
		const view3<double> pf = Component<F>::p_nonlocal_2();
		const view3<double> ps = Component<S>::p_nonlocal_2();
		auto cross = [&](const view3<double> &pb, int B) {
			return shifted(pb, x, A, 1) - shifted(pb, x, A, 0) - shifted(pb, x, A, 1, B, -1) + shifted(pb, x, B, -1);
		};
		return inv_d(A) * inv_d(A) * (shifted(pa, x, A, 1) - 2.0 * shifted(pa, x, A, 0) + shifted(pa, x, A, -1)) +
			   inv_d(A) * inv_d(F) * cross(pf, F) + inv_d(A) * inv_d(S) * cross(ps, S);
	}

	/**
	 * @brief Slot m when it takes the Drude + critical points stretch update, or -1
	 *
	 * A nonlinear material, SHG or THG run, takes the same Drude + critical points update:
	 * its nonlinear polarisation is added to P beforehand, nonlinear::update_polarisation, as in the per-cell path.
	 * The nonlocal model and the electron-density DFT keep the per-cell path.
	 */
	static int stretch_index(int m)
	{
		const bool ok = is(m, Kind::drude_cp) && !(fdtd::materials::model_spec().nonlocal && !fdtd::materials::model_spec().hydrodynamics);
		return ok ? m : -1;
	}

	/**
	 * @brief Drude + critical points ADE at cell k of a row, local model
	 * @verbatim
	      E_new  = c0 curl + c1 (-a0 P_prev + (1 - a1) P) + c2 E_prev + c3 E + sum_n c1 (-b0_n Pn_prev + (1 - b1_n) Pn)
	               + boundary terms
	      P_new  = a0 P_prev + a1 P + a2 (E_prev + 2 E + E_new)
	      Pn_new = b0_n Pn_prev + b1_n Pn + b2_n E_prev + b3_n E + b4_n E_new
	   @endverbatim
	 * The one copy of this arithmetic: the stretch update calls it with no boundary terms, the per-cell
	 * update, materials::models::Drude2CP, inc/materials/models/drude_2cp.h, with the CPML / TF-SF terms.
	 * @param cf        the material's coefficients, materials::models::Drude2CP::Coeffs
	 * @param boundary  boundary(c0, e_new) adds the CPML and TF/SF terms scaled by c0
	 */
	template <class Boundary>
	__attribute__((always_inline)) static void point(const Coeffs &cf, int mm, int k, const stretch::AdeRow &r, double rotdiff, Boundary &&boundary)
	{
		const auto &mp = fdtd::materials::model_parameters();
		const double c0 = cf.c0;
		const double c1 = cf.c1;
		const double c2 = cf.c2;
		const double c3 = cf.c3;
		const double a0 = cf.a0;
		const double a1 = cf.a1;
		const double a2 = cf.a2;
		const int nr_cp = mp.nr_cp[mm];
		double e_new = c0 * rotdiff + c1 * (-a0 * r.p_prev[k] + (1.0 - a1) * r.p[k]) +
					  c2 * r.e_prev[k] +
					  c3 * r.e[k];
		for (int n = 0; n < nr_cp; ++n)
			e_new += c1 * (-cf.b0[n] * r.pc_prev[n][k] + (1.0 - cf.b1[n]) * r.pc[n][k]);
		boundary(c0, e_new);
		const double p_new = a0 * r.p_prev[k] + a1 * r.p[k] + a2 * (r.e_prev[k] + 2.0 * r.e[k] + e_new);
		double pc_new[max_poles];
		for (int n = 0; n < nr_cp; ++n)
		{
			pc_new[n] = cf.b0[n] * r.pc_prev[n][k] + cf.b1[n] * r.pc[n][k] + cf.b2[n] * r.e_prev[k] +
						cf.b3[n] * r.e[k] + cf.b4[n] * e_new;
		}
		r.e_prev[k] = r.e[k];
		r.e[k] = e_new;
		r.p_prev[k] = r.p[k];
		r.p[k] = p_new;
		for (int n = 0; n < nr_cp; ++n)
			r.pc_prev[n][k] = r.pc[n][k], r.pc[n][k] = pc_new[n];
	}

	/**
	 * @brief Drude + critical points ADE over [k0, k1) of row (i, j) of component A
	 * @param boundary  boundary(k, c0, e_new) adds the CPML psi and TF/SF terms of cell k, if any
	 */
	template <int A, class Curl, class Boundary>
	__attribute__((always_inline)) static void stretch(int mm, int i, int j, int k0, int k1, Curl &&curl, Boundary &&boundary)
	{
		using C = Component<A>;
		const Coeffs &cf = coeffs<Drude2CP>(mm);
		const stretch::AdeRow r = stretch::make_ade_row(C::e(), C::e_prev(), C::p(), C::p_prev(), C::pc(), C::pc_prev(), mm, i, j);
		for (int k = k0; k < k1; ++k)
			point(cf, mm, k, r, curl(k), [&](double c, double &e_new) { boundary(k, c, e_new); });
	}

	template <int A, class Boundary>
	static void update(const Cell &c, Boundary &&boundary)
	{
		const auto &mp = fdtd::materials::model_parameters();
		using namespace fdtd::constants; // EP0, MU0, ...
		using fdtd::grid::subdomain;
		using C = Component<A>;
		const view3<double> e = C::e();
		const view3<double> e_prev = C::e_prev();
		auto &pc = C::pc();
		auto &pc_prev = C::pc_prev();
		const int mm = c.mm;
		const int i = c.i;
		const int j = c.j;
		const int k = c.k;
		const int ii = i + 1;
		const int jj = j + 1;
		const int kk = k + 1;
		const int nr_cp = mp.nr_cp[mm];
		const Coeffs &cf = coeffs<Drude2CP>(mm);
		const bool nonlocal = fdtd::materials::model_spec().nonlocal && !fdtd::materials::model_spec().hydrodynamics;
		if (!nonlocal)
		{
			const stretch::AdeRow r = stretch::make_ade_row(e, e_prev, C::p(), C::p_prev(), pc, pc_prev, mm, i, j);
			point(cf, mm, k, r, c.curl, boundary);
			return;
		}
		// nonlocal model: P lives in the padded p_nonlocal, and grad(div P) enters both updates
		const view3<double> p_nl = C::p_nonlocal();
		const view3<double> p_nl_prev = C::p_nonlocal_prev();
		const view3<double> p_nl_2 = C::p_nonlocal_2();
		double Grad_Div = 0;
		double GNOR_Correct = 0;
		double P_CP_new[max_poles];
		const bool inner = c.J2 < fdtd::materials::hydro::nonlocal_barrier_int && (j >= subdomain.j_cpml_min && j <= subdomain.j_cpml_max) &&
						   (i >= subdomain.i_cpml_min && i <= subdomain.i_cpml_max) && (k >= subdomain.k_cpml_min && k <= subdomain.k_cpml_max);
		if (gradient_divergence && inner)
			Grad_Div = grad_div<A>({ii, jj, kk});
		else if (laplacian && inner)
			Grad_Div = fdtd::grid::lattice().inv_d[0] * fdtd::grid::lattice().inv_d[0] * (p_nl_2[ii + 1][jj][kk] - 2.0 * p_nl_2[ii][jj][kk] + p_nl_2[ii - 1][jj][kk]) +
					   fdtd::grid::lattice().inv_d[1] * fdtd::grid::lattice().inv_d[1] * (p_nl_2[ii][jj + 1][kk] - 2.0 * p_nl_2[ii][jj][kk] + p_nl_2[ii][jj - 1][kk]) +
					   fdtd::grid::lattice().inv_d[2] * fdtd::grid::lattice().inv_d[2] * (p_nl_2[ii][jj][kk + 1] - 2.0 * p_nl_2[ii][jj][kk] + p_nl_2[ii][jj][kk - 1]);
		double E_new = cf.c0 * c.curl + cf.c1 * (-cf.a0 * p_nl_prev[i][j][k] + (1.0 - cf.a1) * p_nl[ii][jj][kk]) +
					   cf.c2 * e_prev[i][j][k] +
					   cf.c3 * e[i][j][k] - cf.c0 * fdtd::materials::hydro::non_loc[mm] * Grad_Div / fdtd::materials::model_spec().dt;
		if (fdtd::materials::model_spec().gnor)
		{
			GNOR_Correct = fdtd::materials::hydro::gnor_1[mm] * Grad_Div + fdtd::materials::hydro::gnor_2[mm] * (3.0 * Grad_Div - 4.0 * C::grad_div_n()[i][j][k] + C::grad_div_n_1()[i][j][k]);
			C::grad_div()[i][j][k] = Grad_Div;
			E_new -= cf.c0 * GNOR_Correct / fdtd::materials::model_spec().dt;
		}
		for (int n = 0; n < nr_cp; ++n)
			E_new += cf.c1 * (-cf.b0[n] * pc_prev[mm][n][i][j][k] + (1.0 - cf.b1[n]) * pc[mm][n][i][j][k]);
		boundary(cf.c0, E_new);

		double P_new = cf.a0 * p_nl_prev[i][j][k] +
					   cf.a1 * p_nl[ii][jj][kk] +
					   cf.a2 * (e_prev[i][j][k] + 2.0 * e[i][j][k] + E_new) + fdtd::materials::hydro::non_loc[mm] * Grad_Div;
		if (fdtd::materials::model_spec().gnor)
			P_new += GNOR_Correct;
		for (int n = 0; n < nr_cp; ++n)
		{
			P_CP_new[n] = cf.b0[n] * pc_prev[mm][n][i][j][k] +
						  cf.b1[n] * pc[mm][n][i][j][k] +
						  cf.b2[n] * e_prev[i][j][k] +
						  cf.b3[n] * e[i][j][k] +
						  cf.b4[n] * E_new;
		}
		e_prev[i][j][k] = e[i][j][k];
		e[i][j][k] = E_new;
		p_nl_prev[i][j][k] = p_nl[ii][jj][kk];
		p_nl[ii][jj][kk] = P_new;
		C::p()[mm][i][j][k] = P_new;
		for (int n = 0; n < nr_cp; ++n)
		{
			pc_prev[mm][n][i][j][k] = pc[mm][n][i][j][k];
			pc[mm][n][i][j][k] = P_CP_new[n];
		}
	}
};
} // namespace fdtd::materials::models

#endif
