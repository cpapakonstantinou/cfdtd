/**
 * @file materials/models/drude.h
 * @brief Drude model, codes 150-199
 */
#ifndef FDTD_MATERIALS_MODELS_DRUDE_H
#define FDTD_MATERIALS_MODELS_DRUDE_H

#include "materials/component.h"
#include "materials/models/cell.h"
#include "materials/coefficients.h"
#include "materials/stretch.h"
#include "materials/models/dispersive_e.h"
#include "materials/spec.h"

namespace fdtd::materials::models
{
/**
 * @brief Drude ADE
 * @verbatim
      P_new = (4 P + (beta - 2) P_prev + 2 eps0 alpha E) / (beta + 2)
      E     = ce1 E + ce2 curl - ce2 (P_new - P) / dt
   @endverbatim
 * The CPML psi and TF/SF terms of this model are not added here: the E-face psi pass, boundary::e_face
 * and the plane-wave excitation, planeWave.cc, add them after the scan, scaled by ce_2_d.
 */
struct Drude
{
	/// update coefficients of one material
	struct Coeffs
	{
		double alpha, beta; ///< omega_p^2 dt^2 and gamma dt
		double ce1, ce2;    ///< E update, DispersiveE
	};

	/// coefficients of dispersive material index mm, material code `code`, for time step dt
	static Coeffs compute(int mm, int code, double dt)
	{
		const auto &mp = fdtd::materials::model_parameters();
		using namespace fdtd::constants; // EP0, MU0, ...
		const DispersiveE e = dispersive_e(mm, code, dt);
		return {mp.ome_plasma[mm] * mp.ome_plasma[mm] * dt * dt, mp.gam_plasma[mm] * dt, e.ce1, e.ce2};
	}

	/// the update at cell k of a row, the 3D model and the 1D auxiliary grid; no boundary terms
	static void point(const Coeffs &cf, int k, const stretch::AdeRow &r, double curl)
	{
		const double p_new = (4.0 * r.p[k] + (cf.beta - 2.0) * r.p_prev[k] + 2.0 * constants::EP0 * cf.alpha * r.e[k]) / (cf.beta + 2.0);
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
		const Coeffs &cf = coeffs<Drude>(mm);
		const stretch::AdeRow r{e[i][j], nullptr, C::p()[mm][i][j], C::p_prev()[mm][i][j], {}, {}};
		for (int k = k0; k < k1; ++k)
			point(cf, k, r, curl(k));
	}

	template <int A, class Boundary>
	static void update(const Cell &c, Boundary &&)
	{
		using C = Component<A>;
		const view3<double> e = C::e();
		const int mm = c.mm, i = c.i, j = c.j;
		point(coeffs<Drude>(mm), c.k, stretch::AdeRow{e[i][j], nullptr, C::p()[mm][i][j], C::p_prev()[mm][i][j], {}, {}}, c.curl);
	}
};
} // namespace fdtd::materials::models

#endif
