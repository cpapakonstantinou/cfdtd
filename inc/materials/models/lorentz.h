/**
 * @file materials/models/lorentz.h
 * @brief Lorentz model, model lorentz
 */
#ifndef FDTD_MATERIALS_MODELS_LORENTZ_H
#define FDTD_MATERIALS_MODELS_LORENTZ_H

#include "materials/component.h"
#include "materials/models/cell.h"
#include "materials/coefficients.h"
#include "materials/stretch.h"
#include "materials/models/dispersive_e.h"
#include "materials/spec.h"

namespace fdtd::materials::models
{
/**
 * @brief Lorentz ADE
 * @verbatim
      P_new = ((2 - alpha) P - (1 - beta) P_prev + eps0 (eps_s - eps_inf) alpha E) / (1 + beta)
      E     = ce1 E + ce2 curl - ce2 (P_new - P) / dt
   @endverbatim
 * With hydrodynamics the new P is also stored in pc_sum, padded index.
 * The CPML psi and TF/SF terms of this model are not added here: the E-face psi pass, boundary::e_face
 * and the plane-wave excitation, planeWave.cc, add them after the scan, scaled by ce_2_d.
 */
struct Lorentz
{
	/// update coefficients of one material
	struct Coeffs
	{
		double alpha; ///< omega_p^2 dt^2 and gamma dt
		double beta; ///< omega_p^2 dt^2 and gamma dt
		double ce1; ///< E update, DispersiveE
		double ce2; ///< E update, DispersiveE
	};

	/// coefficients of material slot mm for time step dt
	static Coeffs compute(int mm, double dt)
	{
		const auto &mp = fdtd::materials::model_parameters();
		using namespace fdtd::constants; // EP0, MU0, ...
		const DispersiveE e = dispersive_e(mm, dt);
		return {mp.ome_plasma[mm] * mp.ome_plasma[mm] * dt * dt, mp.gam_plasma[mm] * dt, e.ce1, e.ce2};
	}

	/// the update at cell k of a row, the 3D model and the 1D auxiliary grid; no boundary terms
	static void point(const Coeffs &cf, int mm, int k, const stretch::AdeRow &r, double curl)
	{
		const auto &mp = fdtd::materials::model_parameters();
		using namespace fdtd::constants; // EP0, MU0, ...
		const double p_new = ((2.0 - cf.alpha) * r.p[k] - (1.0 - cf.beta) * r.p_prev[k] + EP0 * (mp.epre_sta[mm] - mp.epre_inf[mm]) * cf.alpha * r.e[k]) / (1.0 + cf.beta);
		r.e[k] *= cf.ce1;
		r.e[k] += cf.ce2 * curl;
		r.e[k] -= cf.ce2 * (p_new - r.p[k]) / fdtd::materials::model_spec().dt;
		r.p_prev[k] = r.p[k];
		r.p[k] = p_new;
	}

	/// cells [k0, k1) of row (i, j): the row is set up once, curl(k) is the cell's curl
	template <int A, class Curl>
	static void stretch(int mm, int i, int j, int k0, int k1, Curl &&curl)
	{
		using C = Component<A>;
		const view3<double> e = C::e();
		const Coeffs &cf = coeffs<Lorentz>(mm);
		const stretch::AdeRow r{e[i][j], nullptr, C::p()[mm][i][j], C::p_prev()[mm][i][j], {}, {}};
		const bool hydro = fdtd::materials::model_spec().hydrodynamics;
		for (int k = k0; k < k1; ++k)
		{
			point(cf, mm, k, r, curl(k));
			if (hydro)
				C::pc_sum()[i + 1][j + 1][k + 1] = r.p[k];
		}
	}

	template <int A, class Boundary>
	static void update(const Cell &c, Boundary &&)
	{
		using C = Component<A>;
		const view3<double> e = C::e();
		const int mm = c.mm;
		const int i = c.i;
		const int j = c.j;
		const int k = c.k;
		const stretch::AdeRow r{e[i][j], nullptr, C::p()[mm][i][j], C::p_prev()[mm][i][j], {}, {}};
		point(coeffs<Lorentz>(mm), mm, k, r, c.curl);
		if (fdtd::materials::model_spec().hydrodynamics)
			C::pc_sum()[i + 1][j + 1][k + 1] = r.p[k];
	}
};
} // namespace fdtd::materials::models

#endif
