/**
 * @file materials/coefficients.cc
 * @brief The update coefficients of the materials (compute_ade_coefficients)
 */

#include "materials/materials.h"
#include "materials/models.h"
#include "grid/lattice.h"
using fdtd::constants::E0;
using fdtd::constants::EP0;
using fdtd::constants::MU0;

namespace
{
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace

/**
 * @brief The update coefficients of every material for the run's time step
 *
 * Dielectrics: ce_1, ce_2 and the refractive index. Dispersive materials: each model computes its own
 * coefficients, `Model::compute`, inc/materials/models/, into its table, materials/coefficients.h;
 * the nonlocal and electron-gas constants of the hydro module are set here. Then the H coefficients.
 */
void fdtd::materials::compute_ade_coefficients()
{
	auto &mp = fdtd::materials::model_parameters();
	const double dt = fdtd::materials::model_spec().dt;

	for (int m = 0; m < slot_count(); ++m)
	{
		const int n = m;
		switch (kind(m))
		{
		case Kind::dielectric:
			mp.ce_1[m] = (double(2.) * EP0 * mp.epre[m] - dt * mp.sig_m[m]) / (double(2.) * EP0 * mp.epre[m] + dt * mp.sig_m[m]);
			mp.ce_2[m] = double(2.) * dt / (double(2.) * EP0 * mp.epre[m] + dt * mp.sig_m[m]);
			mp.n_refr_index[m] = sqrt(double(mp.epre[m]));
			break;
		case Kind::drude_cp:
		{
			coeff_table<models::Drude2CP>[n] = models::Drude2CP::compute(n, dt);
			const double beta = mp.gam_plasma[n] * dt;
			const bool nonlocal = fdtd::materials::model_spec().nonlocal;
			fdtd::materials::hydro::non_loc[n] = nonlocal ? 2 * dt * dt * fdtd::materials::hydro::beta_nl[n] * fdtd::materials::hydro::beta_nl[n] / (2 + beta) : 0;
			fdtd::materials::hydro::gnor_1[n] = nonlocal ? 2 * dt * dt * fdtd::materials::hydro::diff_coeff[n] * mp.gam_plasma[n] / (2 + beta) : 0;
			fdtd::materials::hydro::gnor_2[n] = nonlocal ? 2 * dt * dt * fdtd::materials::hydro::diff_coeff[n] / (2 * dt * (2 + beta)) : 0;
			break;
		}
		case Kind::drude:
			coeff_table<models::Drude>[n] = models::Drude::compute(n, dt);
			break;
		case Kind::lorentz:
			coeff_table<models::Lorentz>[n] = models::Lorentz::compute(n, dt);
			break;
		case Kind::modified_lorentz:
			coeff_table<models::ModifiedLorentz>[n] = models::ModifiedLorentz::compute(n, dt);
			break;
		case Kind::hydrodynamic:
		case Kind::drude_lorentz:
			// the electron gas of the hydrodynamic and Drude-Lorentz models
			fdtd::materials::hydro::powerterm[n] = pow(fdtd::materials::hydro::equilibrium_electron_density[n], 2.0 / 3.0);
			fdtd::materials::hydro::inverse_edens[n] = 1.0 / fdtd::materials::hydro::equilibrium_electron_density[n];
			fdtd::materials::hydro::inverse_chargedens[n] = 1.0 / (E0 * fdtd::materials::hydro::equilibrium_electron_density[n]);
			if (kind(m) == Kind::hydrodynamic)
				coeff_table<models::Hydrodynamic>[n] = models::Hydrodynamic::compute(n, dt);
			else
				coeff_table<models::DrudeLorentz>[n] = models::DrudeLorentz::compute(n, dt);
			break;
		case Kind::magnetic_lorentz:
			coeff_table<models::MagneticLorentz>[n] = models::MagneticLorentz::compute(n, dt);
			break;
		}
	}

}
