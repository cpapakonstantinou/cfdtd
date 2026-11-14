/**
 * @file sources/spec.h
 * @brief What the sources need from the run: time step, background media, the structure centre, the
 *        band, which source is on and its parameters, the output folders
 *
 * Set by the application once the run is configured, end of app::initialise_run; the sources read
 * it instead of the run configuration, docs/dev/Roadmap.md §2.5.
 */
#ifndef FDTD_SOURCES_SPEC_H
#define FDTD_SOURCES_SPEC_H

#include <ctime>
#include <string>

namespace fdtd::sources
{
struct Spec
{
	// time
	double dt = 0;               ///< time step [s]
	double coeff_h = 0;          ///< H update coefficient dt / mu0, the 1D grid
	int one_dim_steps = 0;       ///< steps of the standalone 1D simulations
	int one_dim_cells = 0;       ///< cells of the 1D grid
	// background and structure
	int medium[3] = {};          ///< first, middle and second, last, medium, material codes
	double middle_layer = 0;     ///< middle layer thickness [m]
	int centre[3] = {};          ///< the structure's centre cell
	bool grating = false;        ///< the periodic round_hybrid grating case
	bool hydrodynamics = false;  ///< hydrodynamic electron gas, the dipole source reads the density
	bool harmonic = false;       ///< separate harmonic fields (*_nl)
	// band
	bool shg = false, thg = false;   ///< harmonic generation: the upper half of the band is the harmonic
	int nr_frequencies = 0;          ///< analysed frequencies
	double fundamental_frequency = 0, highest_frequency = 0; ///< [Hz]
	double start_cw = 0;             ///< CW turn-on [periods]
	// which source and its parameters
	bool tfsf = false, mode_source = false, dipole = false, gauss_beam = false, superposition = false;
	bool oblique_with_medium = false; ///< oblique plane waves through the layered background
	double circ_time_shift = 0;       ///< delay between the two components of a circular wave [s]
	int signal_type = 0;              ///< pulse kind code
	bool cw_pulse = false, continuous_wave = false; ///< the CW signals
	int theta_intervals = 0, phi_intervals = 0;     ///< Gaussian beam: plane-wave samples
	double theta_max_degrees = 0, phi_max_degrees = 0;
	double waist_factor = 0, focal_factor = 0;      ///< Gaussian beam
	// outputs
	int plot_xy = 0, plot_xz = 0, plot_yz = 0;      ///< frame periods of the plot planes, 0 off
	std::string reports_dir, input_dir;
	clock_t log_start = 0;                          ///< program start, for log timestamps
};
extern Spec the_spec; ///< set by set_spec
inline const Spec &spec() { return the_spec; }
inline void set_spec(const Spec &s) { the_spec = s; }
} // namespace fdtd::sources

#endif
