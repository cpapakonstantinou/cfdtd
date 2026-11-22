/**
 * @file monitors/spec.h
 * @brief What the monitors need from the run: time step and step count, the band, the analyses
 *        asked for, the structure and background, the source, the output folders
 *
 * Set by the application once the step count is final, Driver::configure_run, after the pulse
 * transform; the monitors read it instead of the run configuration, docs/dev/Roadmap.md §2.5.
 */
#ifndef FDTD_MONITORS_SPEC_H
#define FDTD_MONITORS_SPEC_H

#include <string>

namespace fdtd::monitors
{
struct Spec
{
	// time and band
	double dt = 0;                   ///< time step [s]
	int steps = 0;                   ///< the run's step count
	int nr_frequencies = 0;          ///< analysed frequencies
	bool shg = false; ///< the upper half of the band is the harmonic
	bool thg = false; ///< the upper half of the band is the harmonic
	bool harmonic = false;           ///< separate harmonic fields (*_nl)
	bool hydrodynamics = false;      ///< hydrodynamic electron gas, Q coefficients add the hydro terms
	// analyses
	bool abs_sigma = false;
	bool art = false;
	bool sphere_multipoles = false;
	bool dft_plot = false;
	bool save_output = false;
	int nr_qsca[3] = {};             ///< scattering flux box inset [cells]
	// structure and background
	int shape = 0;                   ///< [geometry] shape code
	bool draw = false;               ///< a shape is configured
	bool graphene_optim = false;
	bool grating = false;
	int centre[3] = {};              ///< structure centre cell
	double antenna_shift = 0; ///< [m]
	double thickness = 0; ///< [m]
	int antenna_material = 0;
	int shell_material = 0;
	int gap_material = -1;           ///< the anapole's gap fill, -1 for none
	int medium_first = 0;
	int medium_second = 0;
	// source
	double e_peak = 0;
	double fundamental_frequency = 0;
	int polarization_deg = 0;
	// output
	long global_size[3] = {};        ///< output block size along x, y, z, inner cells + 1; 0 on CPML-only ranks
	std::string fields_dir;
	std::string frames_dir;
	std::string reports_dir;
};
extern Spec the_spec; ///< set by set_spec
inline const Spec &spec() { return the_spec; }
inline void set_spec(const Spec &s) { the_spec = s; }
} // namespace fdtd::monitors

#endif
