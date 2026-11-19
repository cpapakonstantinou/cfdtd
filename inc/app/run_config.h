/**
 * @file app/run_config.h
 * @brief The run configuration: values fixed once the run is configured, the same on every rank
 *
 * Filled from SimConfig by `apply_sim_config` and completed by the configure phases of the driver, `apply_sim_config`, `app::setup_domain`, `app::initialise_run`,
 * then read-only: `run_config()` gives const access, `edit_run_config()` throws once the driver has sealed it.
 * Grouped by topic so a module can be handed the part it needs. Names are the legacy global names for now.
 */
#ifndef FDTD_APP_RUN_CONFIG_H
#define FDTD_APP_RUN_CONFIG_H

#include "materials/definition.h"
#include <string>
#include <vector>
#include "geometry/shapes.h"
#include "sources/kinds.h"

namespace fdtd::app
{
/// Values fixed for the whole run, see the file comment.
struct RunConfig
{
	/// Time stepping
	struct Time
	{
		double dt;                                                 ///< Time step in s.
		double coeff_h;                                            ///< H update coefficient DT / MU0, the 1/d factors are applied in the update.
		double cfl_coeff;                                          ///< Courant factor: DT = cfl_coeff * ddd / C0.
		int one_dim_nr_max_time_steps;                             ///< Number of steps of the standalone 1D simulation.
	} time;

	/// Global grid, the same on every rank
	struct Grid
	{
		double dx;                                                 ///< Cell size along x in m.
		double dy;                                                 ///< Cell size along y in m.
		double dz;                                                 ///< Cell size along z in m.
		int num_cpml_cells_x;                                      ///< CPML thickness on each x face in cells.
		int num_cpml_cells_y;                                      ///< CPML thickness on each y face in cells.
		int num_cpml_cells_z;                                      ///< CPML thickness on each z face in cells.
		bool pbc_x;                                                ///< Periodic boundaries along x.
		bool pbc_z;                                                ///< Periodic boundaries along z.
		int nr_cells_uni;                                          ///< 1D grid: number of cells (pw_factor * nr_cells_tot_y).
		int nr_cells_tot_x;                                        ///< Global cells along x, CPML included, as configured, adjusted for periodic plane-wave runs.
		int nr_inflate_x;                                          ///< TF/SF box inset from the inner domain along x in cells.
		int nr_cells_tot_y;                                        ///< Global cells along y, CPML included, as configured, adjusted for periodic plane-wave runs.
		int nr_inflate_y;                                          ///< TF/SF box inset from the inner domain along y in cells.
		int nr_cells_tot_z;                                        ///< Global cells along z, CPML included, as configured, adjusted for periodic plane-wave runs.
		int nr_inflate_z;                                          ///< TF/SF box inset from the inner domain along z in cells.
	} grid;

	/// Source kind and signal
	struct Source
	{
		bool dipole_source;                                        ///< True when the excitation is a set of point dipoles, source.type = dipole.
		bool mode_source;                                          ///< Excitation is a waveguide mode read from file.
		bool tfsf;                                                 ///< Total-field/scattered-field plane-wave source is active.
		bool continuous_wave;                                      ///< Source signal: continuous wave.
		bool cw_pulse;                                             ///< Source signal: CW with a smooth turn-on.
		double start_cw;                                           ///< Turn-on time of the CW source in s.
		double fundamental_frequency;                              ///< Carrier frequency of the fundamental in Hz.
		double highest_frequency;                                  ///< Highest frequency of the spectrum in Hz.
		int pw_factor;                                             ///< Length of the 1D auxiliary grid as a multiple of the 3D y extent.
		bool plane_wave_superposition = false;                     ///< Excitation is a superposition of oblique plane waves.
		bool gauss_beam = false;                                   ///< Excitation is a Gaussian beam.
		bool oblique_incidence_with_medium;                        ///< Oblique plane waves cross a layered background.
		bool enforce_one_dim_simulation;                           ///< Run the 1D simulations of the superposed waves even without a layered background.
		bool circ_lh;                                              ///< Left-handed circular polarisation.
		bool circ_rh;                                              ///< Right-handed circular polarisation.
		double circ_time_shift;                                    ///< Time shift between the two linear components of a circular source in s.
		double e_peak;                                             ///< Peak amplitude of the source field in V/m.
		double chirp;                                              ///< Linear chirp of the Gaussian signal.
		double chirp_2;                                            ///< Quadratic chirp of the Gaussian signal.
		int cw_periods;                                            ///< Length of the CW pulse [periods].
		double pulse_width;                                         ///< Source pulse width in s.
		double center_freq;                                        ///< Centre wavelength of the pulse in m, despite the name: f_0 = C0 / center_freq.
		int kind;                                                  ///< Source kind code, `source::Kind`; `setup_domain` maps the periodic plane wave to code 1.
		int polarization_deg;                                      ///< Polarisation angle of the incident wave in deg: 0 x, 90 z.
		bool allow_truncated_signal;                               ///< Run even when the time steps do not cover the signal.
		int signal_type;                                           ///< Signal kind code, `pulse::Kind`; 4 for the mode source.
		int theta_nr_intervals;                                    ///< Gaussian beam: polar angle samples.
		int phi_nr_intervals;                                      ///< Gaussian beam: azimuth samples.
		double theta_max_degrees;                                  ///< Gaussian beam: largest polar angle in deg.
		double phi_max_degrees;                                    ///< Gaussian beam: largest azimuth in deg.
		double beam_w_0_coeff;                                     ///< Gaussian beam: waist factor.
		double focal_length_coeff;                                 ///< Gaussian beam: focal factor.
		int tfsf_y_shrink;                                         ///< TF/SF box shrink along y in cells; set by `setup_domain` for periodic plane waves.
	} source;

	/// Materials and electron models
	struct Materials
	{
		std::vector<fdtd::materials::Definition> definitions; ///< the run's materials, slot s = definitions[s], slot 0 vacuum
		bool magnetic_lorentz;                                     ///< The antenna uses the magnetic Lorentz model.
		int first_medium;                                          ///< Material slot of the background on the source side, low y.
		int middle_medium;                                         ///< Material slot of the middle layer between first_medium and second_medium.
		int second_medium;                                         ///< Material slot of the background on the far side, high y.
		int antenna_material;                                      ///< Material slot of the antenna.
		bool nonlocal;                                             ///< Electrons are nonlocal, electrons != "local".
		bool gnor;                                                 ///< Use the GNOR, diffusive, nonlocal model.
		bool hydrodynamics;                                        ///< The hydrodynamic electron model is active.
		double fermi_v;                                            ///< [electrons] fermi_velocity, the factor on every material's Fermi velocity.
	} materials;

	/// Nonlinear optics
	struct Nonlinear
	{
		bool pnl_approach;                                         ///< Nonlinear update mode: true = add P_NL to the polarisation, false = add it to the E update.
		bool second_harmonic;                                                  ///< Second-harmonic generation is active.
		bool third_harmonic;                                                  ///< Third-harmonic generation is active.
		bool seperate_linear_nonlinear;                            ///< Keep the harmonic fields in separate arrays (*_nl).
		int shg;                                                   ///< Second-harmonic generation: 0 off, 1 on.
		int thg;                                                   ///< Third-harmonic generation: 0 off, 1 on.
	} nonlinear;

	/// Two-temperature model
	struct Thermal
	{
		bool pump_probe;                                           ///< Run a pump simulation followed by a probe simulation.
		bool pump;                                                 ///< True when the two-temperature model is driven by a pump pulse.
		bool probe;                                                ///< This run is the probe of a pump-probe pair.
		bool save_output_temperature;                              ///< Write temperature frames during the run.
		bool load_temperature;                                     ///< Load a static temperature map that sets the Drude damping.
		bool two_temperature_model;                                ///< The two-temperature model is active.
	} thermal;

	/// Geometry: the shape, `antenna_type`, tested with `is(shape::...)`, and its dimensions
	struct Geometry
	{
		/// The configured shape is `k`, `shape::Kind`.
		bool is(int k) const { return antenna_type == k; }
		bool graphene_optim;                                       ///< Average |E| in the SiO2 under a graphene layer, optimisation output.
		bool sta_ir_casing;                                        ///< Wrap the antenna in an IR casing.
		bool nano_antenna;                                         ///< An antenna shape is present, as opposed to a bare background.
		bool grating;                                              ///< Shape flag: grating.
		int center_position_x;                                     ///< Global x index of the shape centre.
		int center_position_y;                                     ///< Global y index of the shape centre.
		int center_position_z;                                     ///< Global z index of the shape centre.
		double width;                                              ///< Antenna width in m.
		double thickness;                                          ///< Antenna thickness along y in m.
		double length;                                             ///< Antenna arm length in m.
		double length_2;                                            ///< Second arm length for asymmetric shapes in m.
		double middle_layer;                                       ///< Thickness of the middle layer in m.
		double rotation_deg;                                       ///< In-plane rotation of the shapes that take one in deg.
		double antenna_shift;                                      ///< Shift of the antenna along y relative to the interface in m.
		int shell_material;                                         ///< Material code of the shell of a core-shell sphere.
		double length_1_l;                                         ///< Composite antennas: length along z in m.
		double length_2_b;                                         ///< Composite antennas: length along x in m.
		int antenna_type;                                          ///< Shape code, `shape::Kind`.
	} geometry;

	/// Spectra and analysis outputs
	struct Analysis
	{
		int nr_an_frequencies;                                     ///< Number of analysed frequencies.
		bool dft_plot;                                             ///< Write the DFT plane fields.
		double sweep_lambda_long;                                  ///< Longest wavelength of the spectrum in m.
		double sweep_lambda_short;                                 ///< Shortest wavelength of the spectrum in m.
		bool abs_sigma;                                            ///< Compute absorption from sigma |E|^2.
		bool art;                                                  ///< Compute absorption, reflection and transmission.
		bool multipoles;                                           ///< Enables the multipole decomposition of the scattered field.
		int nr_qsca_x;                                             ///< Scattering flux box inset along x in cells.
		int nr_qtot_x;                                             ///< Total, absorption, flux box inset along x in cells.
		int nr_qsca_y;                                             ///< Scattering flux box inset along y in cells.
		int nr_qtot_y;                                             ///< Total, absorption, flux box inset along y in cells.
		int nr_qsca_z;                                             ///< Scattering flux box inset along z in cells.
		int nr_qtot_z;                                             ///< Total, absorption, flux box inset along z in cells.
		int dft_plane_mask;                                        ///< Display DFT planes as a bit mask: 1 yz, 2 xz, 4 xy, 8 with H.
	} analysis;

	/// Paths and output switches
	struct Output
	{
		std::string base_path_input;                               ///< Directory holding the input data files, materials, chi, modes.
		std::string base_path_output;                              ///< Root output directory of the run.
		bool save_output;                                          ///< Master switch for writing field/DFT output files.
		int plot_yz;                                               ///< yz animation frame period, steps; 0 off.
		int plot_xz;                                               ///< xz animation frame period, steps; 0 off.
		int plot_xy;                                               ///< xy animation frame period, steps; 0 off.
		int time_plot_3d;                                          ///< 3D field frame period, steps; 0 off.
		int dft_phase_norm;                                        ///< Normalise the DFT by the source spectrum phase (0/1).
		std::string reports_dir;                                    ///< <output>/reports: text tables, Q, ART, signal, point files, 1D profiles, and logs.
		std::string fields_dir;                                     ///< <output>/fields: frequency-domain data, DFT planes, boxes, volumes, near-to-far faces.
		std::string frames_dir;                                     ///< <output>/frames: time-domain frames and probe traces.
		std::string materials_dir;                                  ///< <output>/materials: material index maps, magnetisation, the geometry report.
		std::string state_dir;                                      ///< <output>/state: application state dumps, temperatures, restart.
	} output;
};

/// The run configuration, const. Inline: hot loops read it like the globals it replaces.
extern RunConfig the_run_config;
inline const RunConfig &run_config() { return the_run_config; }

/// Write access for the configure phases; throws std::logic_error after seal_run_config().
RunConfig &edit_run_config();
/// Called by the driver when configuration is complete.
void seal_run_config();
void reopen_run_config();            ///< `fdtd memory` only: the configure phase runs once per virtual rank
} // namespace fdtd::app

#endif
