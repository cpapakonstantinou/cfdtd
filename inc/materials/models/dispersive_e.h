/**
 * @file materials/models/dispersive_e.h
 * @brief The E-update coefficients of the single-pole dispersive models, Drude, Lorentz, magnetic Lorentz
 */
#ifndef FDTD_MATERIALS_MODELS_DISPERSIVE_E_H
#define FDTD_MATERIALS_MODELS_DISPERSIVE_E_H

#include "materials/state.h"
#include "core/constants.h"

namespace fdtd::materials::models
{
/// E = ce1 E + ce2 curl - ce2 (P_new - P) / dt, from eps_inf and the conductivity
struct DispersiveE
{
	double ce1, ce2;
};

/// the coefficients of dispersive material index mm, material code `code`, for time step dt
inline DispersiveE dispersive_e(int mm, int code, double dt)
{
	const auto &mp = fdtd::materials::model_parameters();
	using namespace fdtd::constants; // EP0, MU0, ...
	const double uu1 = mp.sig_m[code] * dt / double(2.);
	const double uu2 = EP0 * mp.epre_inf[mm] + uu1;
	return {(EP0 * mp.epre_inf[mm] - uu1) / uu2, dt / uu2};
}
} // namespace fdtd::materials::models

#endif
