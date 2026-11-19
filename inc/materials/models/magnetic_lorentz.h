/**
 * @file materials/models/magnetic_lorentz.h
 * @brief Magnetic, magneto-optic, Lorentz model, model magnetic_lorentz
 */
#ifndef FDTD_MATERIALS_MODELS_MAGNETIC_LORENTZ_H
#define FDTD_MATERIALS_MODELS_MAGNETIC_LORENTZ_H

#include "materials/component.h"
#include "materials/models/cell.h"
#include "materials/coefficients.h"
#include <cmath>
#include "materials/models/dispersive_e.h"
#include "materials/spec.h"

namespace fdtd::materials::models
{
/**
 * @brief Magnetic Lorentz ADE
 * @verbatim
      cross = M_{B1} dE_{B2} - M_{B2} dE_{B1},   M averaged to the E_A point, dE = E_scratch - E_prev2
              Ex: My dEz - Mz dEy,  Ey: Mz dEx - Mx dEz,  Ez: Mx dEy - My dEx
      P_new = c0 P + c1 P_prev + c2 P_prev2 + c3 P_prev3 + b0 E + b1 E_prev + b2 E_prev2 + b3 cross
      E     = ce1 E + ce2 curl - ce2 (P_new - P) / dt
   @endverbatim
 * With hydrodynamics the new P is also stored in pc_sum, padded index. Neither this update nor the
 * E-face pass adds CPML psi or TF/SF terms for this model, so it is only usable inside the TF/SF box.
 */
struct MagneticLorentz
{
	/// update coefficients of one material: the P recursion, c0..c3 on P, P_prev, P_prev_2, P_prev_3; b0..b3 of pole 0, and the E update
	struct Coeffs : Poles
	{
		double ce1; ///< E update, DispersiveE
		double ce2; ///< E update, DispersiveE
	};

	/// coefficients of material slot mm for time step dt
	static Coeffs compute(int mm, double dt)
	{
		const auto &mp = fdtd::materials::model_parameters();
		using namespace fdtd::constants; // EP0, MU0, ...
		const int n = mm;
		double cp2;
		Coeffs c{};
		const DispersiveE e = dispersive_e(mm, dt);
		c.ce1 = e.ce1;
		c.ce2 = e.ce2;
		cp2 = (1.0 / (dt*dt*dt*dt) + mp.gam_plasma[n]/(dt*dt*dt));
		c.c0 = (4.0/(dt*dt*dt*dt) + 2.0*mp.gam_plasma[n]/(dt*dt*dt) - (2.0*mp.ome_plasma[n]*mp.ome_plasma[n]+mp.gam_plasma[n]*mp.gam_plasma[n])/(dt*dt) - mp.ome_plasma[n]*mp.ome_plasma[n]*mp.gam_plasma[n]/dt)/cp2;
		c.c1 = (-mp.ome_plasma[n]*mp.ome_plasma[n]*mp.ome_plasma[n]*mp.ome_plasma[n] - 6.0/(dt*dt*dt*dt) + 2.0*(2.0*mp.ome_plasma[n]*mp.ome_plasma[n]+mp.gam_plasma[n]*mp.gam_plasma[n])/(dt*dt))/cp2;
		c.c2 = (4.0/(dt*dt*dt*dt) - 2.0*mp.gam_plasma[n]/(dt*dt*dt) - (2.0*mp.ome_plasma[n]*mp.ome_plasma[n]+mp.gam_plasma[n]*mp.gam_plasma[n])/(dt*dt) + mp.ome_plasma[n]*mp.ome_plasma[n]*mp.gam_plasma[n]/dt)/cp2;
		c.c3 = (- 1.0 / (dt*dt*dt*dt) + mp.gam_plasma[n]/(dt*dt*dt))/cp2;

		c.b0[0] = mp.ome_plasma[n]*mp.ome_plasma[n]*EP0*(mp.epre_sta[n] - mp.epre_inf[n])*(1.0/(dt*dt) + mp.gam_plasma[n]/(2.0*dt))/cp2;
		c.b1[0] = mp.ome_plasma[n]*mp.ome_plasma[n]*EP0*(mp.epre_sta[n] - mp.epre_inf[n])*(mp.ome_plasma[n]*mp.ome_plasma[n] - 2.0/(dt*dt))/cp2;
		c.b2[0] = mp.ome_plasma[n]*mp.ome_plasma[n]*EP0*(mp.epre_sta[n] - mp.epre_inf[n])*(1.0/(dt*dt) - mp.gam_plasma[n]/(2.0*dt))/cp2;
		c.b3[0] = -EP0*mp.a_3[n]*mp.ome_plasma[n]/(2.0*dt)/cp2;
		return c;
	}

	template <int A, class Boundary>
	static void update(const Cell &c, Boundary &&)
	{
		using namespace fdtd::constants; // EP0, MU0, ...
		using C = Component<A>;
		using C1 = Component<C::B1>;
		using C2 = Component<C::B2>;
		const view3<double> e = C::e();
		const view3<double> e_prev = C::e_prev();
		const view3<double> e_prev_2 = C::e_prev_2();
		auto &p = C::p();
		auto &p_prev = C::p_prev();
		auto &p_prev_2 = C::p_prev_2();
		auto &p_prev_3 = C::p_prev_3();
		const int mm = c.mm;
		const int i = c.i;
		const int j = c.j;
		const int k = c.k;
		const Coeffs &cf = coeffs<MagneticLorentz>(mm);
		const int x[3] = {i, j, k};
		auto m_at_e = [&](const view3<double> &m) { return 0.5 * m[i][j][k] + 0.5 * shifted(m, x, A, 1); };
		auto de = [&](const view3<double> &scratch, const view3<double> &prev_2) { return scratch[i][j][k] - prev_2[i][j][k]; };
		const double M1 = m_at_e(C1::magnetization());
		const double M2 = m_at_e(C2::magnetization());
		const double cross = M1 * de(C2::e_scratch(), C2::e_prev_2()) - M2 * de(C1::e_scratch(), C1::e_prev_2());

		const double p_new = cf.c0 * p[mm][i][j][k] + cf.c1 * p_prev[mm][i][j][k] + cf.c2 * p_prev_2[mm][i][j][k] + cf.c3 * p_prev_3[mm][i][j][k] +
							 cf.b0[0] * e[i][j][k] + cf.b1[0] * e_prev[i][j][k] + cf.b2[0] * e_prev_2[i][j][k] + cf.b3[0] * cross;
		e[i][j][k] *= cf.ce1;
		e[i][j][k] += cf.ce2 * c.curl;
		e[i][j][k] -= cf.ce2 * (p_new - p[mm][i][j][k]) / fdtd::materials::model_spec().dt;
		p_prev_3[mm][i][j][k] = p_prev_2[mm][i][j][k];
		p_prev_2[mm][i][j][k] = p_prev[mm][i][j][k];
		p_prev[mm][i][j][k] = p[mm][i][j][k];
		p[mm][i][j][k] = p_new;
		if (fdtd::materials::model_spec().hydrodynamics)
			C::pc_sum()[i + 1][j + 1][k + 1] = p_new;
	}
};
} // namespace fdtd::materials::models

#endif
