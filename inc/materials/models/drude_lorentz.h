/**
 * @file materials/models/drude_lorentz.h
 * @brief Drude + one Lorentz pole, codes 350-399; "nonlinear" refers to the electron gas
 *
 * What this register is, checked 2026-11-13 against the code and its history:
 * - The Drude part is the free-electron current J. It is the nonlinear part: with the hydrodynamic
 *   model on, J comes from the hydrodynamic solver, src/hydro, and with the two-temperature model
 *   wp and gamma follow the electron temperature. Otherwise it is a plain linear Drude term.
 * - The Lorentz part is one linear pole, P'' + G P' + W^2 P = eps0 A^2 E, stored in the critical-point
 *   arrays, big_omega_cp = W, big_gamma_cp = G, a_coeff_cp = A = W sqrt(delta eps); the phase phi_cp
 *   is not used. It is the left side of Varin's nonlinear Lorentz equation with the linear drive
 *   only, C. Varin et al., Comput. Phys. Commun. 222, 70, 2018; Opt. Mater. Express 9, 771, 2019:
 *       P'' + G P' + W^2 P = W^2 eps0 (chi1 E + chi2 E^2 + chi3 E^3)
 *   The chi2 / chi3 drive is not implemented. The legacy code set chi_nl = 2.16e-18 m^2/V^2 and
 *   computed NONLINLorentz1 = beta / cp1 and NONLINLorentz2 = (mu0 e / m_e) / cp1, but no update
 *   ever used them; the SHG/THG tensor path, nl_code, src/nonlinear, is not enabled for these codes
 *   either, tensor.cc accepts dispersive slots below ModLorentz only.
 * - Adding the drive is the natural SHG / THG path for these materials: chi2 E^2 gives the second
 *   harmonic, chi3 E^3 the third harmonic and the Kerr term, both with the pole's dispersion. To stay
 *   explicit the drive must use E at step n only, Varin's leapfrog form; the E_new term of the
 *   averaged linear drive below cannot carry it. Roadmap: Phase 2, nonlinear review.
 * - Slots: 350 is ITO, pole present only for middle_thickness < 50 nm; thicker films are pure
 *   Drude, 351 a custom material from the read_in_* keys.
 */
#ifndef FDTD_MATERIALS_MODELS_DRUDE_LORENTZ_H
#define FDTD_MATERIALS_MODELS_DRUDE_LORENTZ_H

#include "materials/thermal/thermal.h"
#include "materials/component.h"
#include "materials/models/cell.h"
#include "materials/coefficients.h"
#include <cmath>
#include "materials/thermal/state.h"
#include "materials/spec.h"

namespace fdtd::materials::models
{
/**
 * @brief Nonlinear Drude-Lorentz ADE
 * @verbatim
      Drude current J (P_new = -J):
        hydrodynamics:              J from the hydrodynamic solver, padded index
        two-temperature / loaded T: J = (1/dt - g/2)/(1/dt + g/2) J + eps0 wp(Te)^2 / (1/dt + g/2) E,  g = g(Te, Tl)
        otherwise:                  J = d1 J + d2 E
      E_new   = c0 (curl + P_new) + c2 E_prev + c3 E + c1 (-b0 P_L,prev + (1 - b1) P_L) + CPML/TF-SF terms * c0
      P_L,new = b0 P_L,prev + b1 P_L + b2 E_prev + b3 E + b4 E_new, Lorentz pole 0
   @endverbatim
 */
struct DrudeLorentz
{
	/// update coefficients of one material: one pole and the electron gas (nl1, nl2)
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
		Coeffs c{};
		NL1 = 1.0 / dt + 0.5 * mp.gam_plasma[n];
		c.nl1 = (1.0 / dt - 0.5 * mp.gam_plasma[n]) / NL1;

		if (fdtd::materials::model_spec().two_temperature || fdtd::materials::model_spec().load_temperature)
		{
			c.nl2 = EP0 / NL1;
		}
		else
		{
			c.nl2 = mp.ome_plasma[n] * mp.ome_plasma[n] * EP0 / NL1;
		}

		for (int s = 0; s < 1; ++s)
		{
			cp1 = (1.0 / (dt * dt) + mp.big_gamma_cp[n][s] / (2.0 * dt) + mp.big_omega_cp[n][s] * mp.big_omega_cp[n][s] / 4.0);
			c.b0[s] = -(1.0 / (dt * dt) - mp.big_gamma_cp[n][s] / (2.0 * dt) + mp.big_omega_cp[n][s] * mp.big_omega_cp[n][s] / 4.0) / cp1;
			c.b1[s] = -(-2.0 / (dt * dt) + mp.big_omega_cp[n][s] * mp.big_omega_cp[n][s] / 2.0) / cp1;
			c.b2[s] = EP0 * mp.a_coeff_cp[n][s] * mp.a_coeff_cp[n][s] / (4.0 * cp1);
			c.b3[s] = EP0 * mp.a_coeff_cp[n][s] * mp.a_coeff_cp[n][s] / (2.0 * cp1);
			c.b4[s] = EP0 * mp.a_coeff_cp[n][s] * mp.a_coeff_cp[n][s] / (4.0 * cp1);
		}

		cp2 = (1.0 / dt + c.b4[0] / (dt * EP0 * mp.epre_inf[n]));
		c.c0 = 1.0 / (cp2 * EP0 * mp.epre_inf[n]);
		c.c1 = 1.0 / (dt * cp2 * EP0 * mp.epre_inf[n]);
		c.c2 = -c.b2[0] / (cp2 * dt * EP0 * mp.epre_inf[n]);
		c.c3 = (1.0 / dt - c.b3[0] / (dt * EP0 * mp.epre_inf[n])) / cp2;
		return c;
	}

	template <int A, class Boundary>
	static void update(const Cell &c, Boundary &&boundary)
	{
		using namespace fdtd::constants; // EP0, MU0, ...
		using C = Component<A>;
		const view3<double> e = C::e();
		const view3<double> e_prev = C::e_prev();
		auto &jc = C::j();
		auto &pc = C::pc();
		auto &pc_prev = C::pc_prev();
		const int mm = c.mm;
		const int i = c.i;
		const int j = c.j;
		const int k = c.k;
		const Coeffs &cf = coeffs<DrudeLorentz>(mm);
		const int ii = i + 1;
		const int jj = j + 1;
		const int kk = k + 1;
		double P_new;
		if (!fdtd::materials::model_spec().hydrodynamics)
		{
			if (fdtd::materials::model_spec().two_temperature || fdtd::materials::model_spec().load_temperature)
			{
				const double gamma = fdtd::materials::thermal::damping_constant(fdtd::materials::thermal::electron_temp[i][j][k], fdtd::materials::thermal::lattice_temp[i][j][k]);
				jc[mm][ii][jj][kk] = ((1.0 / fdtd::materials::model_spec().dt - 0.5 * gamma) / (1.0 / fdtd::materials::model_spec().dt + 0.5 * gamma)) * jc[mm][ii][jj][kk] + (EP0 * fdtd::materials::thermal::plasma_frequency(fdtd::materials::thermal::electron_temp[i][j][k]) * fdtd::materials::thermal::plasma_frequency(fdtd::materials::thermal::electron_temp[i][j][k]) / (1.0 / fdtd::materials::model_spec().dt + 0.5 * gamma)) * e[i][j][k];
				P_new = -jc[mm][ii][jj][kk];
			}
			else
			{
				jc[mm][i][j][k] = cf.nl1 * jc[mm][i][j][k] + (cf.nl2) * e[i][j][k];
				P_new = -jc[mm][i][j][k];
			}
		}
		else
			P_new = -jc[mm][ii][jj][kk];
		double E_new = cf.c0 * (c.curl + P_new) + cf.c2 * e_prev[i][j][k] + cf.c3 * e[i][j][k];
		E_new += cf.c1 * (-cf.b0[0] * pc_prev[mm][0][i][j][k] + (1.0 - cf.b1[0]) * pc[mm][0][i][j][k]);
		boundary(cf.c0, E_new);
		const double P_CP_new = cf.b0[0] * pc_prev[mm][0][i][j][k] +
								cf.b1[0] * pc[mm][0][i][j][k] +
								cf.b2[0] * e_prev[i][j][k] +
								cf.b3[0] * e[i][j][k] +
								cf.b4[0] * E_new;
		pc_prev[mm][0][i][j][k] = pc[mm][0][i][j][k];
		pc[mm][0][i][j][k] = P_CP_new;
		e_prev[i][j][k] = e[i][j][k];
		e[i][j][k] = E_new;
	}
};
} // namespace fdtd::materials::models

#endif
