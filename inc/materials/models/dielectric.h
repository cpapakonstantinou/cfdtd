/**
 * @file materials/models/dielectric.h
 * @brief A dispersionless material, model dielectric: E = c1 E + c2 curl, over whole stretches of a row
 */
#ifndef FDTD_MATERIALS_MODELS_DIELECTRIC_H
#define FDTD_MATERIALS_MODELS_DIELECTRIC_H

#include "materials/slots.h"
#include "materials/nonlinear/state.h"
#include "materials/spec.h"

namespace fdtd::materials::models
{
struct Dielectric
{
	/// a plain dielectric, neither dispersive nor nonlinear; seperate_linear_nonlinear and mode_source runs
	/// keep the per-cell path, which also updates the nonlinear field / uses the full material map
	static bool applies(int m)
	{
		return !fdtd::materials::model_spec().cellwise && !fdtd::materials::model_spec().harmonic && !dispersive(m) &&
			   !(fdtd::materials::nonlinear::nl[m] && (fdtd::materials::nonlinear::iso_chi_2[m] || fdtd::materials::nonlinear::iso_chi_3[m]));
	}

	/// e = c1 e + c2 curl(k) over [k0, k1)
	template <class Curl>
	__attribute__((always_inline)) static void stretch(double *e, int k0, int k1, double c1, double c2, Curl &&curl)
	{
		for (int k = k0; k < k1; ++k)
		{
			const double rotdiff = curl(k);
			e[k] *= c1;
			e[k] += c2 * rotdiff;
		}
	}
};
} // namespace fdtd::materials::models

#endif
