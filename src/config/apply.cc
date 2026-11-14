/** \file   config/apply.cc
 *  \brief  Maps the typed SimConfig onto app::RunConfig and the remaining solver globals.
 *
 *  Replaces ReadJSONDataFile() / ReadIniDataFile() and the side files they
 *  pulled in, pphinfoCTW.txt, num_rk.txt, PumpProbeDelay.txt, Probes.txt,
 *  pphREADINMATERIAL.txt, pphnlsini.txt, pphmatini.txt,
 *  gaussian_beam.txt, the chi2/ and chi3/ tensor files.
 *
 *  Unlike the legacy readers this runs on EVERY rank, from the same INI text,
 *  so globals that were only ever set on the root process, update_n,
 *  with_crit_point, ..., are now consistent everywhere.
 **/
#include <cctype>
#include "app/run_config.h"
#include "app/run_state.h"
#include "app/structures.h"
#include "io/binary.h"
#include "sources/one_dim_simulator.h"
#include "io/text.h"
#include "materials/thermal/state.h"
#include "grid/state.h"
#include "materials/state.h"
#include "geometry/state.h"
#include "monitors/state.h"
#include "sources/state.h"
#include "materials/hydro/state.h"
#include "config/apply.h"
#include "app/structures.h"
#include "sources/mode_source.h"
#include "geometry/render.h"
#include "config/exceptions.h"

using fdtd::constants::C0;
using fdtd::constants::PI2;

namespace
{
fdtd::app::RunConfig &run_cfg = fdtd::app::edit_run_config(); ///< configure phase: write access
} // namespace

namespace fdtd::config
{
namespace
{

/// Copy the tensor called `name`, unset = zero, into `dst`, scaled: the first `count` components.
void copy_tensor(int order, const std::string& name, const std::vector<ChiTensor>& user, double* dst, int count, double scale)
{
	const ChiTensor* t = name.empty() ? nullptr : find_chi_tensor(order, name, user);
	if (!name.empty() && !t)
		throw validation_error("unknown chi" + std::to_string(order) + " tensor '" + name + "'");
	for (int i = 0; i < count; ++i)
		dst[i] = t ? double(scale * t->values[i]) : double(0);
}

std::string lower(std::string s)
{
	for (auto& ch : s)
		ch = char(std::tolower(static_cast<unsigned char>(ch)));
	return s;
}

/**
 * Sets the solver's model flags from [electrons], [hydrodynamic], [temperature], [two_temperature]. Reproduces the former "Advanced Model"
 * codes exactly:
 *   0 local | 1 nonlocal | 2 gnor
 *   3..15   hydrodynamic with the pressure / convection / magnetic / update_density /
 *           critical_point / shg_source combinations of the old table
 *   200/201 hydrodynamic + split_nonlinear, convection+magnetic / all terms
 *   202     temperature = two_temperature;   2020 the same with save_temperature
 *   2021    temperature = load;   2022 / 2023 temperature = pump / probe
 *   203     hydrodynamic without pressure/convection/magnetic + two_temperature + save
 */
void apply_model(const SimConfig& c, bool frames_on)
{
	const HydrodynamicSection& m = c.hydrodynamic;
	const std::string electrons = lower(c.electrons.model), temperature = lower(c.temperature.model);
	const bool hydro = electrons == "hydrodynamic";
	const bool thermal = temperature != "none";

	run_cfg.materials.nonlocal = electrons != "local";
	run_cfg.materials.gnor = electrons == "gnor";
	run_cfg.materials.hydrodynamics = hydro;
	if (!run_cfg.materials.nonlocal && !thermal)
		return; // local Maxwell: nothing else changes

	fdtd::geometry::add_free_electron_layer = false;
	if (hydro)
	{
		if (run_cfg.materials.antenna_material < 17 && run_cfg.materials.second_medium < 17 && run_cfg.materials.middle_medium < 17)
			throw validation_error("[electrons] model = hydrodynamic needs a hydrodynamic material (index >= 17)");
		fdtd::materials::hydro::with_pressure = m.pressure;
		fdtd::materials::hydro::with_convection = m.convection;
		fdtd::materials::hydro::with_mag_field = m.magnetic;
		if (m.split_nonlinear)
			run_cfg.nonlinear.seperate_linear_nonlinear = true;
	}
	if (thermal)
	{
		run_cfg.thermal.two_temperature_model = temperature != "load";
		run_cfg.thermal.load_temperature = temperature == "load";
		run_cfg.thermal.save_output_temperature = c.two_temperature.save_frames && frames_on;
		run_cfg.thermal.pump_probe = temperature == "pump" || temperature == "probe";
		run_cfg.thermal.pump = temperature == "pump";
		run_cfg.thermal.probe = temperature == "probe";
	}
}

} // namespace

void apply_sim_config(const SimConfig& c, int /*world_size*/, bool is_root)
{
	auto &mp = fdtd::materials::model_parameters();
	// domain
	run_cfg.grid.num_cpml_cells_x = c.grid.pml_x;
	run_cfg.grid.num_cpml_cells_y = c.grid.pml_y;
	run_cfg.grid.num_cpml_cells_z = c.grid.pml_z;
	run_cfg.grid.nr_cells_tot_x = c.grid.nx + 2 * c.grid.pml_x;
	run_cfg.grid.nr_cells_tot_y = c.grid.ny + 2 * c.grid.pml_y;
	run_cfg.grid.nr_cells_tot_z = c.grid.nz + 2 * c.grid.pml_z;
	fdtd::grid::num_mpi_processes_x = c.parallel.ranks_x;
	fdtd::grid::num_mpi_processes_y = c.parallel.ranks_y;
	fdtd::grid::num_mpi_processes_z = c.parallel.ranks_z;
	run_cfg.grid.dx = double(c.grid.dx);
	run_cfg.grid.dy = double(c.grid.dy);
	run_cfg.grid.dz = double(c.grid.dz);
	fdtd::app::run_state().nr_max_time_steps = c.grid.time_steps;
	run_cfg.geometry.center_position_x = c.grid.origin_x;
	run_cfg.geometry.center_position_y = c.grid.origin_y;
	run_cfg.geometry.center_position_z = c.grid.origin_z;
	run_cfg.grid.pbc_x = c.grid.periodic_x;
	run_cfg.grid.pbc_z = c.grid.periodic_z;
	run_cfg.geometry.sta_ir_casing = c.grid.staircase ? 1 : 0;
	run_cfg.time.cfl_coeff = double(c.grid.courant);

	run_cfg.grid.nr_inflate_x = c.boxes.tfsf_x;
	run_cfg.grid.nr_inflate_y = c.boxes.tfsf_y;
	run_cfg.grid.nr_inflate_z = c.boxes.tfsf_z;
	run_cfg.analysis.nr_qsca_x = c.boxes.scattering_x;
	run_cfg.analysis.nr_qsca_y = c.boxes.scattering_y;
	run_cfg.analysis.nr_qsca_z = c.boxes.scattering_z;
	run_cfg.analysis.nr_qtot_x = c.boxes.absorption_x;
	run_cfg.analysis.nr_qtot_y = c.boxes.absorption_y;
	run_cfg.analysis.nr_qtot_z = c.boxes.absorption_z;

	// load balancing: every rank plans the same split, nothing is broadcast --
	{
		// rank boundaries planned by rank_starts(): never on a TF/SF plane
		std::string note;
		auto chunks = [&](const char* axis, const std::string& raw, int n, int total, int a, std::vector<int>& arr, bool& equal) {
			arr = rank_starts(axis, total, n, raw, tfsf_planes(c, a), &note);
			equal = false;
		};
		chunks("x", c.parallel.chunks_x, c.parallel.ranks_x, c.grid.nx + 2 * c.grid.pml_x, 0, fdtd::grid::arr_decomp_x, fdtd::grid::equal_decomp_x);
		chunks("y", c.parallel.chunks_y, c.parallel.ranks_y, c.grid.ny + 2 * c.grid.pml_y, 1, fdtd::grid::arr_decomp_y, fdtd::grid::equal_decomp_y);
		chunks("z", c.parallel.chunks_z, c.parallel.ranks_z, c.grid.nz + 2 * c.grid.pml_z, 2, fdtd::grid::arr_decomp_z, fdtd::grid::equal_decomp_z);
		if (is_root && !note.empty())
			fdtd::io::log("run", "%s\n", note.c_str());
	}

	fdtd::app::anapole_input = c.anapole;

	// wave superposition, [wave.*]
	fdtd::sources::plane_wave_input.clear();
	for (const auto& w : c.waves)
	{
		PlaneWave p;
		p.Dir_X = w.direction_x; p.Dir_Y = w.direction_y; p.Dir_Z = w.direction_z;
		p.AmpEx = w.e_x; p.AmpEy = w.e_y; p.AmpEz = w.e_z;
		p.AmpHx = w.h_x; p.AmpHy = w.h_y; p.AmpHz = w.h_z;
		p.UseDerivative = c.wave_superposition.derivative;
		fdtd::sources::plane_wave_input.push_back(p);
	}
	fdtd::sources::plane_wave_dispersion = c.wave_superposition.dispersion_correction;
	fdtd::sources::dipole_input.clear();
	for (const auto& d : c.dipoles)
		fdtd::sources::dipole_input.push_back({d.x, d.y, d.z, double(d.e_x), double(d.e_y), double(d.e_z)});

	// source
	// only the section of the selected signal can differ from its defaults (validate_sim_config)
	const std::string signal = lower(c.source.signal);
	const BandLimitedSignal& band = signal == "raised_cosine_derivative" ? c.signal.raised_cosine_derivative
								  : signal == "rectangular" ? c.signal.rectangular : c.signal.raised_cosine;
	GaussianSignal carrier = c.signal.gaussian;
	if (signal == "gaussian_alt")
	{
		carrier.center_wavelength = c.signal.gaussian_alt.center_wavelength;
		carrier.width = c.signal.gaussian_alt.width;
	}
	run_cfg.source.signal_type = pulse_code(c.source.signal);
	run_cfg.source.mode_source = (lower(c.source.type) == "mode");
	fdtd::sources::mode_input = c.mode;
	fdtd::app::blocks = c.blocks;
	run_cfg.source.kind = source_code(c.source.type);
	fdtd::monitors::max_source_frequency = double(band.max_frequency);
	run_cfg.source.allow_truncated_signal = c.source.allow_truncated_signal;
	run_cfg.source.e_peak = double(c.source.peak_field);
	run_cfg.time.one_dim_nr_max_time_steps = c.source.aux_steps;
	run_cfg.source.center_freq = double(carrier.center_wavelength);
	run_cfg.source.pulse_width = double(carrier.width);
	run_cfg.source.chirp = double(carrier.chirp);
	run_cfg.source.chirp_2 = double(carrier.chirp_quadratic);
	run_cfg.source.cw_periods = c.signal.cw_pulse.periods;
	run_cfg.source.start_cw = double(c.signal.continuous_wave.settle_periods);
	run_cfg.source.pw_factor = c.source.aux_length_factor;
	run_cfg.source.enforce_one_dim_simulation = c.wave_superposition.layered_simulation;
	fdtd::sources::n_start = 0;

	// [plane_wave]: E = cos(a) x + sin(a) z with a = polarization_deg
	const std::string axis = lower(c.plane_wave.polarization_axis);
	run_cfg.source.circ_rh = (axis == "rhcp");
	run_cfg.source.circ_lh = (axis == "lhcp");
	run_cfg.source.circ_time_shift = 0;
	if (axis == "x")
		run_cfg.source.polarization_deg = c.plane_wave.polarization_angle;
	else if (axis == "z")
		run_cfg.source.polarization_deg = 90 - c.plane_wave.polarization_angle;
	else
		run_cfg.source.polarization_deg = 45; // equal parts, the z part a quarter period apart, set below

	run_cfg.source.theta_nr_intervals = c.gaussian_beam.theta_steps;
	run_cfg.source.phi_nr_intervals = c.gaussian_beam.phi_steps;
	run_cfg.source.theta_max_degrees = double(c.gaussian_beam.theta_max);
	run_cfg.source.phi_max_degrees = double(c.gaussian_beam.phi_max);
	run_cfg.source.beam_w_0_coeff = double(c.gaussian_beam.waist_factor);
	run_cfg.source.focal_length_coeff = double(c.gaussian_beam.focal_factor);

	// monitors
	// the solver uses exactly the [spectrum] frequencies, computed with its own C0 and PI2
	{
		const auto& s = c.spectrum;
		const int n = s.points, k = harmonic_order(s.harmonic_band);
		fdtd::sources::spectrum_omega.clear();
		for (double w : spectrum_omegas(s, C0, PI2))
			fdtd::sources::spectrum_omega.push_back(double(w));
		run_cfg.analysis.nr_an_frequencies = int(fdtd::sources::spectrum_omega.size());
		run_cfg.analysis.sweep_lambda_short = double(s.lambda_min);
		run_cfg.analysis.sweep_lambda_long = double(n == 1 ? s.lambda_min : s.lambda_max);
		run_cfg.source.fundamental_frequency = n == 1 ? C0 / run_cfg.analysis.sweep_lambda_short : double(0.5) * (C0 / run_cfg.analysis.sweep_lambda_long + C0 / run_cfg.analysis.sweep_lambda_short);
		run_cfg.source.highest_frequency = k * (C0 / run_cfg.analysis.sweep_lambda_short);
		fdtd::sources::fundamental_step_frequency = n > 1 ? (C0 / run_cfg.analysis.sweep_lambda_short - C0 / run_cfg.analysis.sweep_lambda_long) / (n - 1) : 0;
		fdtd::monitors::harmonic_step_frequency = k * fdtd::sources::fundamental_step_frequency;

		// circular polarisation: the z part leads, rhcp, or lags, lhcp, by a quarter period of the carrier
		if (run_cfg.source.circ_rh || run_cfg.source.circ_lh)
		{
			const bool has_own_carrier = signal == "gaussian" || signal == "gaussian_alt";
			const double f = has_own_carrier ? C0 / run_cfg.source.center_freq : run_cfg.source.fundamental_frequency;
			run_cfg.source.circ_time_shift = (run_cfg.source.circ_rh ? 0.25 : -0.25) / f;
		}
	}
	run_cfg.analysis.dft_plane_mask = dft_plane_mask_bits(c.output);
	run_cfg.output.save_output = c.output.volume_dft;
	run_cfg.output.plot_xy = c.output.frames_xy;
	run_cfg.output.plot_xz = c.output.frames_xz;
	run_cfg.output.plot_yz = c.output.frames_yz;
	run_cfg.output.time_plot_3d = c.output.frames_3d;
	run_cfg.output.dft_phase_norm = (c.output.phase_normalize ? 1 : 0);
	run_cfg.analysis.abs_sigma = c.output.sigma_absorption;

	// multipoles
	fdtd::monitors::exact_multipoles = c.multipoles.exact;
	fdtd::monitors::lwa_multipoles = c.multipoles.lwa;
	run_cfg.analysis.multipoles = c.multipoles.sphere_radius > 0;
	fdtd::monitors::sphere_rad_multipole = c.multipoles.sphere_radius;
	fdtd::monitors::ntheta_multipole = c.multipoles.sphere_theta;
	fdtd::monitors::nphi_multipole = c.multipoles.sphere_phi;

	// geometry
	const GeometrySection g = resolved_geometry(c); // the shape's own section written into the pool
	run_cfg.geometry.antenna_type = shape_code(g.shape);
	run_cfg.geometry.nano_antenna = (run_cfg.geometry.antenna_type != fdtd::geometry::shape::none);
	run_cfg.materials.antenna_material = g.material;
	run_cfg.materials.first_medium = g.first_medium;
	run_cfg.materials.middle_medium = g.middle_medium;
	run_cfg.materials.second_medium = g.last_medium;
	run_cfg.geometry.middle_layer = double(g.middle_thickness);
	run_cfg.geometry.antenna_shift = double(g.shift_y);
	run_cfg.geometry.width = double(g.width);
	run_cfg.geometry.length = double(g.length);
	run_cfg.geometry.thickness = double(g.thickness);
	run_cfg.geometry.length_2 = double(g.length_2);
	fdtd::geometry::antenna_gap = double(g.gap);
	run_cfg.geometry.length_1_l = double(g.arm_length);
	run_cfg.geometry.length_2_b = double(g.arm_width);
	run_cfg.geometry.rotation_deg = double(g.rotation);
	fdtd::geometry::corner_radius = double(g.corner_radius);
	run_cfg.geometry.shell_material = g.shell_material;

	// materials, was pphmatini.txt
	auto n = [&c](int i) {
		auto it = c.materials.refractive_index.find(i);
		return it == c.materials.refractive_index.end() ? double(0) : double(it->second);
	};
	run_cfg.materials.n_refr_0 = n(0);
	run_cfg.materials.n_refr_1 = n(1);
	run_cfg.materials.n_refr_2 = n(2);
	run_cfg.materials.n_refr_3 = n(3);
	for (int i = 4; i < 100; ++i)
		mp.n_refr[i] = n(i);
	run_cfg.materials.carrier_density = double(c.materials.carrier_density);
	run_cfg.materials.carrier_density_back = double(c.materials.carrier_density_back);
	run_cfg.geometry.stretch_eps_x = double(c.materials.stretch_x);

	// user-defined dispersive material, was pphREADINMATERIAL.txt
	const auto& rm = c.custom_material;
	run_cfg.materials.read_in_epre_inf = double(rm.eps_inf);
	run_cfg.materials.read_in_ome_plasma = double(rm.plasma_frequency);
	run_cfg.materials.read_in_gam_plasma = double(rm.plasma_damping);
	run_cfg.materials.read_in_big_omega_cp = double(rm.cp_frequency);
	run_cfg.materials.read_in_big_omega_cp_2 = double(rm.cp_frequency_2);
	run_cfg.materials.read_in_big_gamma_cp = double(rm.cp_damping);
	run_cfg.materials.read_in_big_gamma_cp_2 = double(rm.cp_damping_2);
	run_cfg.materials.read_in_a_coeff_cp = double(rm.cp_amplitude);
	run_cfg.materials.read_in_a_coeff_cp_2 = double(rm.cp_amplitude_2);
	run_cfg.materials.read_in_phi_cp = double(rm.cp_phase);
	run_cfg.materials.read_in_phi_cp_2 = double(rm.cp_phase_2);
	run_cfg.materials.read_in_fermi_velocity = double(rm.fermi_velocity);
	run_cfg.materials.read_in_electron_density = double(rm.electron_density);

	// nonlinear, was pphnlsini.txt + chi2/ chi3/ files; now [chi2|chi3.<name>] or presets
	run_cfg.nonlinear.shg = harmonic_mode(c.nonlinear.shg);
	run_cfg.nonlinear.thg = harmonic_mode(c.nonlinear.thg);
	run_cfg.geometry.graphene_optim = c.nonlinear.graphene_optimization ? 1 : 0;
	run_cfg.nonlinear.pnl_approach = (c.nonlinear.update == "polarization") ? 1 : 0;
	{
		const bool shg = process_harmonic(c.nonlinear) == 2;
		const bool thg = process_harmonic(c.nonlinear) == 3;
		if (shg || thg)
		{
			const int order = shg ? 2 : 3;
			const int count = shg ? 18 : 21;
			const double scale = shg ? double(2) : double(1); // chi2 = 2 d_mn
			auto& nl = run_cfg.nonlinear;
			copy_tensor(order, c.nonlinear.chi_structure, c.chi_tensors, shg ? nl.chi2_antenna : nl.chi3_antenna, count, scale);
			copy_tensor(order, c.nonlinear.chi_first, c.chi_tensors, shg ? nl.chi2_first : nl.chi3_first, count, scale);
			copy_tensor(order, c.nonlinear.chi_middle, c.chi_tensors, shg ? nl.chi2_middle : nl.chi3_middle, count, scale);
			copy_tensor(order, c.nonlinear.chi_last, c.chi_tensors, shg ? nl.chi2_second : nl.chi3_second, count, scale);
		}
	}

	// physics models
	run_cfg.materials.diff_c = double(c.electrons.diffusion);
	run_cfg.materials.fermi_v = double(c.electrons.fermi_velocity);
	fdtd::materials::thermal::num_rk = c.two_temperature.rk_stages;
	fdtd::materials::thermal::pump_probe_delay = double(c.pump_probe.delay);

	// defaults that the advanced-model table may switch on
	run_cfg.thermal.two_temperature_model = false;
	run_cfg.materials.nonlocal = false;
	run_cfg.materials.hydrodynamics = false;
	run_cfg.materials.gnor = false;
	run_cfg.thermal.save_output_temperature = false;
	run_cfg.thermal.pump_probe = run_cfg.thermal.pump = run_cfg.thermal.probe = false;
	fdtd::geometry::add_free_electron_layer = false;
	run_cfg.nonlinear.seperate_linear_nonlinear = false;
	apply_model(c, run_cfg.output.plot_yz > 0 || run_cfg.output.plot_xz > 0 || run_cfg.output.plot_xy > 0);

	// hydrodynamic models make the solver compute a 3w band over multi-point spectra:
	// the user must ask for it, the frequency list never changes silently
	if (run_cfg.materials.hydrodynamics && !run_cfg.nonlinear.seperate_linear_nonlinear && c.spectrum.points > 2 &&
		harmonic_order(c.spectrum.harmonic_band) != 3)
		throw validation_error("[electrons] model = hydrodynamic computes a third-harmonic band over several "
									"points: set [spectrum] harmonic_band = third");
	fdtd::materials::gradient_divergence = true;
	fdtd::materials::laplacian = false;

}

} // namespace fdtd::config
