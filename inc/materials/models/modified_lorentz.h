/**
 * @file materials/models/modified_lorentz.h
 * @brief Modified Lorentz model, model modified_lorentz
 */
#ifndef FDTD_MATERIALS_MODELS_MODIFIED_LORENTZ_H
#define FDTD_MATERIALS_MODELS_MODIFIED_LORENTZ_H

#include "materials/component.h"
#include "materials/models/cell.h"
#include "materials/coefficients.h"
#include "materials/stretch.h"
#include "materials/spec.h"

namespace fdtd::materials::models
{
/**
 * @brief Modified Lorentz ADE
 * @verbatim
      E_new   = ce1 E_prev + ce2 E + ce3 (curl - 1/2 sum_{n=0,1} ((1 + a_n) P_n + b_n P_n,prev)) + CPML/TF-SF terms * ce3
      P_n,new = a_n P_n + b_n P_n,prev + (alpha+_n E_new + beta_n E + alpha-_n E_prev) / dt
   @endverbatim
 * The E update takes the first two poles, the models are fitted with two.
 */
struct ModifiedLorentz
{
	/// update coefficients of one material
	struct Coeffs
	{
		double a[max_poles]; ///< pole recursion: on P_n and on P_n,prev
		double b[max_poles]; ///< pole recursion: on P_n and on P_n,prev
		double alpha_p[max_poles]; ///< pole drive: on E_new and on E_prev
		double alpha_m[max_poles]; ///< pole drive: on E_new and on E_prev
		double beta[max_poles];                      ///< pole drive: on E
		double ce1; ///< E update: on E_prev, on E, on the curl
		double ce2; ///< E update: on E_prev, on E, on the curl
		double ce3; ///< E update: on E_prev, on E, on the curl
	};

	/// coefficients of material slot mm for time step dt
	static Coeffs compute(int mm, double dt)
	{
		const auto &mp = fdtd::materials::model_parameters();
		using namespace fdtd::constants; // EP0, MU0, ...
		const int n = mm;
		Coeffs c{};
		for (int s = 0; s < mp.nr_cp[n]; ++s)
		{
			c.a[s] = (double(2.) - (mp.big_omega_cp[n][s] * dt) * (mp.big_omega_cp[n][s] * dt)) / (mp.big_gamma_cp[n][s] * dt + double(1.));
			c.b[s] = (mp.big_gamma_cp[n][s] * dt - double(1.)) / (mp.big_gamma_cp[n][s] * dt + double(1.));
			c.alpha_p[s] = EP0 * mp.a_coeff_cp[n][s] * dt * (+0.5 * dt * mp.big_omega_cp[n][s] * mp.big_omega_cp[n][s] + mp.phi_cp[n][s]) / (mp.big_gamma_cp[n][s] * dt + 1.);
			c.alpha_m[s] = EP0 * mp.a_coeff_cp[n][s] * dt * (-0.5 * dt * mp.big_omega_cp[n][s] * mp.big_omega_cp[n][s] + mp.phi_cp[n][s]) / (mp.big_gamma_cp[n][s] * dt + 1.);
			c.beta[s] = -2. * EP0 * mp.a_coeff_cp[n][s] * dt * mp.phi_cp[n][s] / (mp.big_gamma_cp[n][s] * dt + 1.);
		}
		const double den = 2. * EP0 * mp.epre_inf[n] + mp.sig_m[n] * dt + c.alpha_p[0] + c.alpha_p[1];
		c.ce1 = -(c.alpha_m[0] + c.alpha_m[1]) / den;
		c.ce2 = (2. * EP0 * mp.epre_inf[n] - mp.sig_m[n] * dt - c.beta[0] - c.beta[1]) / den;
		c.ce3 = 2. * dt / den;
		return c;
	}

	/// the update at cell k of a row, the 3D model and the 1D auxiliary grid; boundary(ce3, e_new)
	/// adds the CPML / TF-SF terms
	template <class Boundary>
	static void point(const Coeffs &cf, int nr_cp, int k, const stretch::AdeRow &r, double curl, Boundary &&boundary)
	{
		double P_CP_new[max_poles];
		double E_new = cf.ce1 * r.e_prev[k] + cf.ce2 * r.e[k] +
					   cf.ce3 * (curl - 0.5 * ((1. + cf.a[0]) * r.pc[0][k] + (1. + cf.a[1]) * r.pc[1][k] + cf.b[0] * r.pc_prev[0][k] + cf.b[1] * r.pc_prev[1][k]));
		boundary(cf.ce3, E_new);
		for (int n = 0; n < nr_cp; ++n)
			P_CP_new[n] = cf.a[n] * r.pc[n][k] + cf.b[n] * r.pc_prev[n][k] + (cf.alpha_p[n] * E_new + cf.beta[n] * r.e[k] + cf.alpha_m[n] * r.e_prev[k]) / fdtd::materials::model_spec().dt;
		r.e_prev[k] = r.e[k];
		r.e[k] = E_new;
		for (int n = 0; n < nr_cp; ++n)
		{
			r.pc_prev[n][k] = r.pc[n][k];
			r.pc[n][k] = P_CP_new[n];
		}
	}

	/// the update over [k0, k1) of row (i, j) of component A; boundary(k, ce3, e_new) adds the CPML /
	/// TF-SF terms of cell k, if any
	template <int A, class Curl, class Boundary>
	static void stretch(int mm, int i, int j, int k0, int k1, Curl &&curl, Boundary &&boundary)
	{
		const auto &mp = fdtd::materials::model_parameters();
		using C = Component<A>;
		const view3<double> e = C::e();
		const view3<double> e_prev = C::e_prev();
		auto &pc = C::pc();
		auto &pc_prev = C::pc_prev();
		const int nr_cp = mp.nr_cp[mm];
		stretch::AdeRow r{e[i][j], e_prev[i][j], nullptr, nullptr, {}, {}};
		for (int n = 0; n < std::max(2, nr_cp); ++n) // the E update reads poles 0 and 1
			r.pc[n] = pc[mm][n][i][j], r.pc_prev[n] = pc_prev[mm][n][i][j];
		const Coeffs &cf = coeffs<ModifiedLorentz>(mm);
		for (int k = k0; k < k1; ++k)
			point(cf, nr_cp, k, r, curl(k), [&](double c, double &e_new) { boundary(k, c, e_new); });
	}

	template <int A, class Boundary>
	static void update(const Cell &c, Boundary &&boundary)
	{
		const auto &mp = fdtd::materials::model_parameters();
		using C = Component<A>;
		const view3<double> e = C::e();
		const view3<double> e_prev = C::e_prev();
		auto &pc = C::pc();
		auto &pc_prev = C::pc_prev();
		const int mm = c.mm;
		const int i = c.i;
		const int j = c.j;
		const int nr_cp = mp.nr_cp[mm];
		stretch::AdeRow r{e[i][j], e_prev[i][j], nullptr, nullptr, {}, {}};
		for (int n = 0; n < std::max(2, nr_cp); ++n) // the E update reads poles 0 and 1
			r.pc[n] = pc[mm][n][i][j], r.pc_prev[n] = pc_prev[mm][n][i][j];
		point(coeffs<ModifiedLorentz>(mm), nr_cp, c.k, r, c.curl, boundary);
	}
};
} // namespace fdtd::materials::models

#endif
