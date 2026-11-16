/**
 * @file materials/models/hydrodynamic.h
 * @brief Hydrodynamic free-electron + critical points model, E update, codes 300-349
 *
 * The hydrodynamic current J comes from the hydrodynamic solver, src/hydrodynamics.cc; this is the
 * E-side point update.
 */
#ifndef FDTD_MATERIALS_MODELS_HYDRODYNAMIC_H
#define FDTD_MATERIALS_MODELS_HYDRODYNAMIC_H

#include "materials/component.h"
#include "io/text.h"
#include "materials/models/cell.h"
#include "materials/coefficients.h"
#include <cmath>
#include "materials/spec.h"

namespace fdtd::materials::models
{
/**
 * @brief Hydrodynamic + critical points E update
 * @verbatim
      P_new  = -J, from the hydrodynamic solver
      E_new  = c0 (curl + P_new) + c2 E_prev + c3 E + sum_n c1 (-b0_n Pn_prev + (1 - b1_n) Pn) + CPML/TF-SF terms * c0
      Pn_new = b0_n Pn_prev + b1_n Pn + b2_n E_prev + b3_n E + b4_n E_new,   pc_sum = sum_n Pn_new
   @endverbatim
 * The separate harmonic field (seperate_linear_nonlinear) takes the same update with its own arrays
 * and curl, without the CPML/TF-SF terms.
 */
struct Hydrodynamic
{
	/// update coefficients of one material: the poles and the electron gas (nl1 .. nl5)
	struct Coeffs : Poles, ElectronGas
	{
	};

	/// coefficients of dispersive material index mm, material code `code`, for time step dt
	static Coeffs compute(int mm, int code, double dt)
	{
		const auto &mp = fdtd::materials::model_parameters();
		using namespace fdtd::constants; // EP0, MU0, ...
		const int n = mm;
		const int m = code;
		(void)m;
		double NL1;
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
		NL1 = 1.0 / dt + 0.5 * mp.gam_plasma[n];
		c.nl1 = (1.0 / dt - 0.5 * mp.gam_plasma[n]) / NL1;
		c.nl2 = 1.0 * E0 * E0 * fdtd::materials::hydro::equilibrium_electron_density[n] / (fdtd::materials::hydro::effective_electron_mass[n] * NL1);
		c.nl3 = 1.0 * E0 * E0 / (fdtd::materials::hydro::effective_electron_mass[n] * NL1);
		c.nl4 = 1.0 / NL1;
		c.nl5 = fdtd::materials::model_spec().fermi_velocity * fdtd::materials::model_spec().fermi_velocity * (E0 * HBAR * HBAR * std::pow(3 * PI * PI, 2.0 / 3.0) / (3.0 * fdtd::materials::hydro::effective_electron_mass[n] * fdtd::materials::hydro::effective_electron_mass[n])) / NL1;

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
		}

		cp9 = mp.sig_m[m] * dt / 2 + c.b3[0] + c.b3[1] - EP0 * mp.epre_inf[n];
		cp10 = mp.sig_m[m] * dt / 2 + c.b4[0] + c.b4[1] + EP0 * mp.epre_inf[n];
		c.c0 = dt / cp10;
		c.c1 = 1 / cp10;
		c.c2 = -(c.b2[0] + c.b2[1]) / cp10;
		c.c3 = -cp9 / cp10;
		return c;
	}

	template <int A, class Boundary>
	static void update(const Cell &c, Boundary &&boundary)
	{
		const auto &mp = fdtd::materials::model_parameters();
		using namespace fdtd::constants; // EP0, MU0, ...
		using C = Component<A>;
		const view3<double> e = C::e();
		const view3<double> e_prev = C::e_prev();
		const view3<double> e_nl = C::e_nl();
		const view3<double> e_prev_nl = C::e_prev_nl();
		const view3<double> pc_2 = C::pc_sum();
		const view3<double> pc_2_nl = C::pc_sum_nl();
		auto &jc = C::j();
		auto &jc_nl = C::j_nl();
		auto &pc = C::pc();
		auto &pc_prev = C::pc_prev();
		auto &pc_nl = C::pc_nl();
		auto &pc_prev_nl = C::pc_prev_nl();
		const int mm = c.mm;
		const int i = c.i;
		const int j = c.j;
		const int k = c.k;
		const Coeffs &cf = coeffs<Hydrodynamic>(mm);
		const int ii = i + 1;
		const int jj = j + 1;
		const int kk = k + 1;
		const int nr_cp = mp.nr_cp[mm];
		const bool nl = fdtd::materials::model_spec().harmonic;
		double P_CP_new[NrMaxCP];
		double P_CP_new_NL[NrMaxCP];
		double E_new_NL = 0;
		const double P_new = -jc[mm][ii][jj][kk];
		double E_new = cf.c0 * (c.curl + P_new) + cf.c2 * e_prev[i][j][k] + cf.c3 * e[i][j][k];
		if (nl)
		{
			const double P_new_NL = -jc_nl[mm][ii][jj][kk];
			E_new_NL = cf.c0 * (c.curl_nl + P_new_NL) + cf.c2 * e_prev_nl[i][j][k] + cf.c3 * e_nl[i][j][k];
		}
		for (int n = 0; n < nr_cp; ++n)
		{
			E_new += cf.c1 * (-cf.b0[n] * pc_prev[mm][n][i][j][k] + (1.0 - cf.b1[n]) * pc[mm][n][i][j][k]);
			if (nl)
				E_new_NL += cf.c1 * (-cf.b0[n] * pc_prev_nl[mm][n][i][j][k] + (1.0 - cf.b1[n]) * pc_nl[mm][n][i][j][k]);
		}
		boundary(cf.c0, E_new);
		pc_2[ii][jj][kk] = 0.0;
		if (nl)
			pc_2_nl[ii][jj][kk] = 0.0;
		for (int n = 0; n < nr_cp; ++n)
		{
			P_CP_new[n] = cf.b0[n] * pc_prev[mm][n][i][j][k] +
						  cf.b1[n] * pc[mm][n][i][j][k] +
						  cf.b2[n] * e_prev[i][j][k] +
						  cf.b3[n] * e[i][j][k] +
						  cf.b4[n] * E_new;
			if (nl)
			{
				P_CP_new_NL[n] = cf.b0[n] * pc_prev_nl[mm][n][i][j][k] +
								 cf.b1[n] * pc_nl[mm][n][i][j][k] +
								 cf.b2[n] * e_prev_nl[i][j][k] +
								 cf.b3[n] * e_nl[i][j][k] +
								 cf.b4[n] * E_new_NL;
				pc_2_nl[ii][jj][kk] += P_CP_new_NL[n];
			}
			pc_2[ii][jj][kk] += P_CP_new[n];
		}
		e_prev[i][j][k] = e[i][j][k];
		e[i][j][k] = E_new;
		if (nl)
		{
			e_prev_nl[i][j][k] = e_nl[i][j][k];
			e_nl[i][j][k] = E_new_NL;
		}
		for (int n = 0; n < nr_cp; ++n)
		{
			pc_prev[mm][n][i][j][k] = pc[mm][n][i][j][k];
			pc[mm][n][i][j][k] = P_CP_new[n];
			if (nl)
			{
				pc_prev_nl[mm][n][i][j][k] = pc_nl[mm][n][i][j][k];
				pc_nl[mm][n][i][j][k] = P_CP_new_NL[n];
			}
		}
	}
};
} // namespace fdtd::materials::models

#endif
