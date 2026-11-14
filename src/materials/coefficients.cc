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

	for (int m = 0; m < NRMAXMATERIALS; ++m)
	{
		if (m < DispMateOffsetD)
		{
			mp.ce_1[m] = (double(2.) * EP0 * mp.epre[m] - dt * mp.sig_m[m]) / (double(2.) * EP0 * mp.epre[m] + dt * mp.sig_m[m]);
			mp.ce_2[m] = double(2.) * dt / (double(2.) * EP0 * mp.epre[m] + dt * mp.sig_m[m]);
			mp.n_refr_index[m] = sqrt(double(mp.epre[m]));
			continue;
		}
		const int n = m - DispMateOffsetD;
		if (n < DrudeCP)
		{
			coeff_table<models::Drude2CP>[n] = models::Drude2CP::compute(n, m, dt);
			const double beta = mp.gam_plasma[n] * dt;
			const bool nonlocal = fdtd::materials::model_spec().nonlocal;
			fdtd::materials::hydro::non_loc[n] = nonlocal ? 2 * dt * dt * fdtd::materials::hydro::beta_nl[n] * fdtd::materials::hydro::beta_nl[n] / (2 + beta) : 0;
			fdtd::materials::hydro::gnor_1[n] = nonlocal ? 2 * dt * dt * fdtd::materials::hydro::diff_coeff[n] * mp.gam_plasma[n] / (2 + beta) : 0;
			fdtd::materials::hydro::gnor_2[n] = nonlocal ? 2 * dt * dt * fdtd::materials::hydro::diff_coeff[n] / (2 * dt * (2 + beta)) : 0;
		}
		else if (n < Drude)
			coeff_table<models::Drude>[n] = models::Drude::compute(n, m, dt);
		else if (n < Lorentz)
			coeff_table<models::Lorentz>[n] = models::Lorentz::compute(n, m, dt);
		else if (n < ModLorentz)
			coeff_table<models::ModifiedLorentz>[n] = models::ModifiedLorentz::compute(n, m, dt);
		else if (n < NonlinearDrudeLorentz)
		{
			// the electron gas of the hydrodynamic and nonlinear Drude-Lorentz models
			fdtd::materials::hydro::powerterm[n] = pow(fdtd::materials::hydro::equilibrium_electron_density[n], 2.0 / 3.0);
			fdtd::materials::hydro::inverse_edens[n] = 1.0 / fdtd::materials::hydro::equilibrium_electron_density[n];
			fdtd::materials::hydro::inverse_chargedens[n] = 1.0 / (E0 * fdtd::materials::hydro::equilibrium_electron_density[n]);
			if (n < NonLinFreeElectron)
				coeff_table<models::Hydrodynamic>[n] = models::Hydrodynamic::compute(n, m, dt);
			else
				coeff_table<models::DrudeLorentz>[n] = models::DrudeLorentz::compute(n, m, dt);
		}
		else if (n < MagneticLorentz)
			coeff_table<models::MagneticLorentz>[n] = models::MagneticLorentz::compute(n, m, dt);
	}

}
