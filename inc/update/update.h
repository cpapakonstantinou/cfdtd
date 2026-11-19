/**
 * @file update/update.h
 * @brief The field updates of one time step
 *
 * @verbatim
      update_h   H from the curl of E on every rank cell, CPML psi on the faces, then the auxiliary 1D grid
      update_e   E per stretch with the material model of the stretch, inc/materials/models, the CPML and
                 TF/SF terms, the nonlinear terms, then the auxiliary 1D grid
   @endverbatim
 * The sources inject after each, sources/source.h; the halo exchange follows in the time loop.
 */
#ifndef FDTD_UPDATE_UPDATE_H
#define FDTD_UPDATE_UPDATE_H

namespace fdtd::update
{
/// What the updates of this run include; set by the application once the run is configured
/// before the first step, read by update_h / update_e instead of the run configuration.
struct Spec
{
	bool tfsf = false;              ///< plane-wave TF/SF corrections
	bool tfsf_xpol = false; ///< the incident wave has an x / z E component, the rows with TF/SF terms
	bool tfsf_zpol = false; ///< the incident wave has an x / z E component, the rows with TF/SF terms
	bool shg = false; ///< second / third harmonic terms
	bool thg = false; ///< second / third harmonic terms
	bool pnl_polarisation = false;  ///< P_NL added to the polarisation, else to the E update
	bool magnetic_lorentz = false;  ///< the magnetic Lorentz model is present
	bool harmonic = false;          ///< separate harmonic fields (*_nl)
	bool cellwise_materials = false; ///< every cell's own material in the CPML, mode source
	int medium_near_y = 0; ///< background material slots below / above the domain in y
	int medium_far_y = 0; ///< background material slots below / above the domain in y
	double coeff_h = 0;             ///< H update coefficient DT / MU0
};
const Spec &spec();
void set_spec(const Spec &s);

void update_h();
void update_e();
} // namespace fdtd::update

#endif
