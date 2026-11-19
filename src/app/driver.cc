/** \file   app/driver.cc
 *  \brief  Phases of a 3D-FDTD run, configuration, decomposition, per-run
 *          preparation, geometry dump, post-processing, output.
 *
 *  The order of calls into the legacy solver is the one of the former
 *  monolithic main(); what changed is that every phase is a method, every
 *  phase reports success collectively, and hardcoded paths / constants come
 *  from the configuration. The time loop lives in time_loop.cc.
 **/
#include "materials/materials.h"
#include "grid/state.h"
#include "materials/thermal/thermal.h"
#include "materials/hydro/hydro.h"
#include "app/setup.h"
#include "io/text.h"
#include "geometry/builders.h"
#include "sources/signal.h"
#include "driver_context.h"
#include "io/h5io.h"
#include "grid/halo_exchange.h"
#include "monitors/custom_boxes.h"
#include "monitors/field_trackers.h"
#include "monitors/time_probes.h"
#include "monitors/volume_dft.h"
#include "sources/waves.h"
#include "config/apply.h"
#include "sources/mode_source.h"
#include "config/errors.h"
#include "materials/stretch.h"
#include "boundary/cpml.h"
#include "core/limits.h"
#include "grid/lattice.h"
#include "app/planes.h"
#include "app/structures.h"
#include "monitors/spec.h"
#include "boundary/incident_field.h"
#include "sources/tfsf.h"
#include "update/update.h"

#include <sys/resource.h>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include "sources/aux_grid.h"
#include "app/run_state.h"
#include "sources/state.h"
#include "materials/thermal/state.h"
#include "monitors/state.h"
#include "materials/nonlinear/state.h"
#include "app/run_config.h"
using fdtd::constants::C0;
using fdtd::constants::PI2;
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;
using fdtd::materials::max_slots;
using fdtd::materials::magnetization_x;
using fdtd::materials::magnetization_y;
using fdtd::materials::magnetization_z;
using fdtd::materials::matind_x;
using fdtd::materials::matind_y;
using fdtd::materials::matind_z;
using fdtd::materials::background_y;

namespace
{
const fdtd::app::RunConfig &run_cfg = fdtd::app::run_config(); ///< the run configuration
} // namespace

namespace fs = std::filesystem;

namespace fdtd
{

// Helpers

namespace
{

/// The output tree below one root: reports/ text tables and logs, fields/ frequency-domain data,
/// frames/ time-domain frames and probe traces, materials/ material index maps, state/ application state dumps.
void assign_output_paths(const std::string& root)
{
	fdtd::app::RunConfig &run_cfg = fdtd::app::edit_run_config(); // configure phase
	run_cfg.output.base_path_output = root;
	run_cfg.output.reports_dir = root + "/reports";
	run_cfg.output.fields_dir = root + "/fields";
	run_cfg.output.frames_dir = root + "/frames";
	run_cfg.output.materials_dir = root + "/materials";
	run_cfg.output.state_dir = root + "/state";
}

std::vector<std::string> output_directories(const std::string& root)
{
	return {root, run_cfg.output.reports_dir, run_cfg.output.fields_dir, run_cfg.output.frames_dir, run_cfg.output.materials_dir, run_cfg.output.state_dir};
}

/// Root reads the file and the catalog files it needs once; every rank gets the same text, no N-way
/// file-system hammering.
bool broadcast_file(const std::string& path, const std::vector<std::string>& overrides, int rank, std::string& text, std::string& err)
{
	long long len = -1;
	if (rank == ROOTPROCESS)
	{
		std::ifstream in(path, std::ios::binary);
		if (in)
		{
			std::ostringstream ss;
			ss << in.rdbuf();
			try
			{
				text = fdtd::config::with_catalog(ss.str(), path, overrides);
				len = static_cast<long long>(text.size());
			}
			catch (const std::exception& e)
			{
				err = e.what();
				len = -2;
			}
		}
	}
	MPI_Bcast(&len, 1, MPI_LONG_LONG, ROOTPROCESS, MPI_COMM_WORLD);
	if (len < 0)
	{
		if (len == -1 || rank != ROOTPROCESS)
			err = len == -1 ? "cannot read configuration file '" + path + "'" : "the materials of the configuration could not be read";
		return false;
	}
	text.resize(static_cast<size_t>(len));
	if (len > 0)
		MPI_Bcast(&text[0], static_cast<int>(len), MPI_CHAR, ROOTPROCESS, MPI_COMM_WORLD);
	return true;
}

/// Peak resident memory of this process [MiB].
double peak_rss_mib()
{
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	return double(ru.ru_maxrss) / 1024.0; // Linux: KiB
}

/// The material index of every E component: the three animation planes, and the volume when asked
void write_material_index(fdtd::grid::OutputLayout &layout, bool volume)
{
	const std::string base = run_cfg.output.materials_dir + "/MaterialIndex";
	int *start = layout.start;
	int *global = layout.global;
	int *local = layout.local;
	int *shift = layout.shift;
	if (volume)
		save_as_binary_3d(base + ".h5", matind_x, matind_y, matind_z, fdtd::grid::output_comm, start, global, local, shift, false, true, false);
	if (subdomain.plots_xy())
		save_as_binary_2d(base + "_XY.h5", matind_x, matind_y, matind_z, fdtd::grid::output_comm_plot_xy, subdomain.k_animation, start, global, local, shift, 0, 1, false, true, false);
	if (subdomain.plots_yz())
		save_as_binary_2d(base + "_YZ.h5", matind_x, matind_y, matind_z, fdtd::grid::output_comm_plot_yz, subdomain.i_animation, start, global, local, shift, 1, 2, false, true, false);
	if (subdomain.plots_xz())
		save_as_binary_2d(base + "_XZ.h5", matind_x, matind_y, matind_z, fdtd::grid::output_comm_plot_xz, subdomain.j_animation, start, global, local, shift, 0, 2, false, true, false);
}
} // namespace

// Construction

Driver::Driver(const CliOptions& opt, int rank, int size)
: _ctx(std::make_unique<Context>()), _opt(opt), _rank(rank), _size(size)
{
}

Driver::~Driver() = default;

bool Driver::agree(bool local_ok, const std::string& what)
{
	int mine = local_ok ? 1 : 0;
	int all = 0;
	int root_ok = mine;
	MPI_Allreduce(&mine, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	MPI_Bcast(&root_ok, 1, MPI_INT, ROOTPROCESS, MPI_COMM_WORLD);
	// a failure on rank 0 is reported once by report_error(); failures only on other
	// ranks, e.g. a file one rank cannot read, are printed by those ranks
	if (!local_ok && root_ok)
		std::fprintf(stderr, "[rank %d] error: %s\n", _rank, what.c_str());
	if (!all && _error.empty())
		_error = what;
	_ok = _ok && (all == 1);
	return _ok;
}

void Driver::report_error()
{
	if (_rank == ROOTPROCESS)
	{
		if (_error.empty())
			std::fprintf(stderr, "\n3D-FDTD %s aborted: see messages above\n", command_name(_opt.command));
		else
			std::fprintf(stderr, "\n3D-FDTD %s aborted: %s\n", command_name(_opt.command), _error.c_str());
	}
}

// Phase: configuration

void Driver::load_configuration()
{
	fdtd::app::RunConfig &run_cfg = fdtd::app::edit_run_config(); // configure phase

	std::string text;
	std::string err;
	if (!agree(broadcast_file(_opt.config_file, _opt.overrides, _rank, text, err), err))
		return;

	// every rank parses the same text -> identical configuration everywhere
	bool good = true;
	try
	{
		_cfg = config::load_sim_config(config::build_config(text, _opt.config_file, _opt.overrides, _opt.input_dir, _opt.output_dir));
		config::validate_sim_config(_cfg);
		const int needed = _cfg.parallel.ranks_x * _cfg.parallel.ranks_y * _cfg.parallel.ranks_z;
		if (_opt.command == Command::Memory)
			_size = fdtd::grid::world_size = needed; // one process plays every rank in turn
		if (needed != _size)
			throw fdtd::config::validation_error("[parallel] ranks_x*ranks_y*ranks_z = " + std::to_string(needed) + " but the job runs on " +
										std::to_string(_size) + " MPI processes");
		const int nfreq = int(config::spectrum_omegas(_cfg.spectrum).size());
		if (nfreq > fdtd::max_frequencies)
			throw fdtd::config::validation_error("[spectrum]: " + std::to_string(nfreq) + " DFT frequencies exceed the limit of " +
										std::to_string(fdtd::max_frequencies));

		run_cfg.output.base_path_input = _cfg.paths.input_dir;
		assign_output_paths(_cfg.paths.output_dir);
		config::apply_sim_config(_cfg, _size, _rank == ROOTPROCESS);
	}
	catch (const std::bad_alloc&)
	{
		good = false;
		err = "out of memory while allocating the configuration structures";
	}
	catch (const std::exception& e)
	{
		good = false;
		err = e.what();
	}
	if (!agree(good, err))
		return;

	_ctx->movie_base = run_cfg.output.frames_dir + "/Efield_";

	// output tree + resolved configuration, root only, result shared with everybody
	good = true;
	if (_rank == ROOTPROCESS && _opt.command != Command::Memory)
	{
		try
		{
			for (const auto& d : output_directories(_cfg.paths.output_dir))
				fs::create_directories(d);
			std::ofstream(_cfg.paths.output_dir + "/resolved_config.ini")
				<< "; resolved configuration of this run (" << command_name(_opt.command) << ", "
				<< _size << " MPI processes, source: " << _opt.config_file << ")\n"
				<< config::to_ini(_cfg);
		}
		catch (const std::exception& e)
		{
			good = false;
			err = std::string("cannot create the output tree: ") + e.what();
		}
		fdtd::io::log("run", "3D-FDTD %s | config %s | output %s | %d MPI processes (%d x %d x %d)\n", command_name(_opt.command),
					_opt.config_file.c_str(), _cfg.paths.output_dir.c_str(), _size, _cfg.parallel.ranks_x, _cfg.parallel.ranks_y, _cfg.parallel.ranks_z);
	}
	if (!agree(good, err))
		return;

}

// `fdtd memory`: play every rank of the layout in turn, up to the allocations before the geometry

bool Driver::estimate_memory(std::vector<std::array<unsigned long long, core::memory::group_count>>& per_rank)
{
	load_configuration();
	if (!_ok)
		return false;
	core::memory::enter_lazy_mode(); // arrays cost memory once written
	fdtd::io::log_muted = true;
	per_rank.assign(_size, {});
	for (int r = 0; r < _size && _ok; ++r)
	{
		fdtd::grid::world_rank = r;
		fdtd::app::reopen_run_config();
		const auto before = core::memory::live;
		const std::size_t held_before = core::memory::held;
		decompose();
		if (_ok)
			initialize_ranks();
		if (_ok)
			configure_run();
		if (_ok)
			allocate_sweep_memory();
		if (_ok)
			build_geometry(); // the dispersive state of structure materials follows their cells
		for (int g = 0; g < core::memory::group_count; ++g)
			per_rank[r][g] = core::memory::live[g] - before[g];
		fdtd::app::release_run_memory();
		monitors::time_probes().release();
		monitors::custom_boxes().release();
		if (_ok && core::memory::held != held_before) // an array outlived the rank: the next rank's count would be off
		{
			_ok = false;
			_error = "internal error: " + std::to_string(core::memory::held - held_before) + " bytes of arrays were not released after rank " +
					 std::to_string(r) + "; use the geometry report of `fdtd mesh` and report this";
			for (const core::memory::Charge* c : core::memory::Charge::outstanding())
				if (c->bytes() > 4096)
					_error += "\n  " + std::string(core::memory::GROUP_NAMES[c->group()]) + ": " + std::to_string(c->bytes()) + " bytes";
		}
	}
	fdtd::grid::world_rank = ROOTPROCESS;
	fdtd::io::log_muted = false;
	return _ok;
}

// Phase: decomposition

void Driver::decompose()
{
	// every rank computes the same decomposition from the same configuration
	fdtd::app::setup_domain();
	{
		const auto &g = fdtd::app::run_config().grid; // final after setup_domain, periodic plane-wave CPML
		fdtd::grid::Lattice lat;
		lat.cells[0] = g.nr_cells_tot_x, lat.cells[1] = g.nr_cells_tot_y, lat.cells[2] = g.nr_cells_tot_z;
		lat.cpml[0] = g.num_cpml_cells_x, lat.cpml[1] = g.num_cpml_cells_y, lat.cpml[2] = g.num_cpml_cells_z;
		lat.tfsf_inset[0] = g.nr_inflate_x, lat.tfsf_inset[1] = g.nr_inflate_y, lat.tfsf_inset[2] = g.nr_inflate_z;
		lat.periodic_x = g.pbc_x, lat.periodic_z = g.pbc_z;
		lat.d[0] = g.dx, lat.d[1] = g.dy, lat.d[2] = g.dz;
		fdtd::grid::set_lattice(lat);
	}
	std::string why;
	const bool good = !fdtd::grid::set_lattice_data(fdtd::app::global_planes(fdtd::grid::lattice()), why);
	if (good)
		fdtd::grid::print_decomposition();
	if (!agree(good, "domain decomposition failed (" + why + ")"))
		return;

}

// Phase: rank initialization

void Driver::initialize_ranks()
{
	fdtd::app::RunConfig &run_cfg = fdtd::app::edit_run_config(); // configure phase
	fdtd::grid::axis_indices(fdtd::grid::world_rank, subdomain.position);

	if (!agree(fdtd::app::initialise_run() && fdtd::sources::build_waves(_ctx->waves, _ctx->frame), "per-rank initialization failed"))
		return;

	auto& c = *_ctx;
	c.layout = grid::output_layout();
	h5io::set_xdmf(_cfg.output.xdmf, run_cfg.grid.dx, run_cfg.grid.dy, run_cfg.grid.dz);

	// the 1D auxiliary simulation is only needed to actually inject the source
	if (run_cfg.source.oblique_incidence_with_medium && !mesh_only())
	{
		if (fdtd::grid::world_rank == ROOTPROCESS)
			fdtd::io::log("run", "Running 1D Simulation\n");
		c.ods.init(c.waves, c.frame);
	}

	fdtd::app::seal_run_config(); // the run configuration is read-only from here on
}

// Phase: prepare the run

void Driver::prepare_run()
{
	_dummy_run = false;
	configure_run();
	if (!_ok)
		return;
	if (!allocate_sweep_memory())
		return;
	build_geometry();
	if (!_ok || mesh_only())
		return;
	setup_sources_and_monitors();
}

/// The pulse spectrum and the per-run flags.
void Driver::configure_run()
{
	std::string why;
	if (!mesh_only() && !run_cfg.source.allow_truncated_signal && !agree(fdtd::sources::signal_fits(fdtd::app::run_state().nr_max_time_steps, why), why))
		return;

	if (_opt.command != Command::Memory) // the signal's spectrum and its files: not needed to size the arrays
		fdtd::app::run_state().nr_max_time_steps = fdtd::sources::transform_pulse(fdtd::app::run_state().nr_max_time_steps);
	{
		monitors::Spec sp; // the step count is final here
		const auto &g = run_cfg.geometry;
		const auto &a = run_cfg.analysis;
		sp.dt = run_cfg.time.dt, sp.steps = fdtd::app::run_state().nr_max_time_steps;
		sp.nr_frequencies = a.nr_an_frequencies;
		sp.shg = run_cfg.nonlinear.second_harmonic, sp.thg = run_cfg.nonlinear.third_harmonic, sp.harmonic = run_cfg.nonlinear.seperate_linear_nonlinear;
		sp.hydrodynamics = run_cfg.materials.hydrodynamics;
		sp.abs_sigma = a.abs_sigma, sp.art = a.art, sp.sphere_multipoles = a.multipoles, sp.dft_plot = a.dft_plot;
		sp.save_output = run_cfg.output.save_output;
		sp.nr_qsca[0] = a.nr_qsca_x, sp.nr_qsca[1] = a.nr_qsca_y, sp.nr_qsca[2] = a.nr_qsca_z;
		sp.shape = g.antenna_type, sp.draw = g.nano_antenna, sp.graphene_optim = g.graphene_optim, sp.grating = g.grating;
		sp.centre[0] = g.center_position_x, sp.centre[1] = g.center_position_y, sp.centre[2] = g.center_position_z;
		sp.antenna_shift = g.antenna_shift, sp.thickness = g.thickness;
		sp.antenna_material = run_cfg.materials.antenna_material, sp.shell_material = g.shell_material, sp.gap_material = 0;
		sp.medium_first = run_cfg.materials.first_medium, sp.medium_second = run_cfg.materials.second_medium;
		sp.e_peak = run_cfg.source.e_peak, sp.fundamental_frequency = run_cfg.source.fundamental_frequency;
		sp.polarization_deg = run_cfg.source.polarization_deg;
		for (int k = 0; k < 3; ++k) // inner cells + 1, 0 on CPML-only ranks
			sp.global_size[k] = subdomain.cpml_only ? 0L : long(fdtd::grid::lattice().inner(k) + 1);
		sp.fields_dir = run_cfg.output.fields_dir, sp.frames_dir = run_cfg.output.frames_dir, sp.reports_dir = run_cfg.output.reports_dir;
		monitors::set_spec(sp);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (run_cfg.analysis.multipoles)
		monitors::field_trackers().setup(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z);
}

bool Driver::allocate_sweep_memory()
{

	// geometry read-in determines which dispersive materials, hence fields, exist

	if (!agree(fdtd::app::allocate_fields(), "cannot allocate the field arrays"))
		return false;
	{
		const int cells[3] = {run_cfg.grid.num_cpml_cells_x, run_cfg.grid.num_cpml_cells_y, run_cfg.grid.num_cpml_cells_z};
		const double du[3] = {run_cfg.grid.dx, run_cfg.grid.dy, run_cfg.grid.dz};
		boundary::cpml().set_profiles(cells, du, run_cfg.time.dt);
	}
	const core::memory::Scope cpml_group(core::memory::cpml);
	const bool cpml_near[3] = {subdomain.nea_cpml_x(), subdomain.nea_cpml_y(), subdomain.nea_cpml_z()};
	const bool cpml_far[3] = {subdomain.far_cpml_x(), subdomain.far_cpml_y(), subdomain.far_cpml_z()};
	if (!agree(boundary::cpml().allocate(subdomain, cpml_near, cpml_far, run_cfg.nonlinear.seperate_linear_nonlinear), "cannot allocate the CPML arrays"))
		return false;
	boundary::set_incident_e_terms(&fdtd::sources::tfsf_e_terms); // the TF/SF surface's incident field
	boundary::cpml().background = {run_cfg.materials.first_medium, run_cfg.materials.second_medium, run_cfg.source.mode_source};
	{
		update::Spec u;
		u.tfsf = run_cfg.source.tfsf;
		u.tfsf_xpol = fdtd::sources::xpol, u.tfsf_zpol = fdtd::sources::zpol;
		u.shg = run_cfg.nonlinear.second_harmonic, u.thg = run_cfg.nonlinear.third_harmonic;
		u.pnl_polarisation = run_cfg.nonlinear.pnl_approach;
		u.magnetic_lorentz = run_cfg.materials.magnetic_lorentz;
		u.harmonic = run_cfg.nonlinear.seperate_linear_nonlinear;
		u.cellwise_materials = run_cfg.source.mode_source;
		u.medium_near_y = run_cfg.materials.first_medium, u.medium_far_y = run_cfg.materials.second_medium;
		u.coeff_h = run_cfg.time.coeff_h;
		update::set_spec(u);
	}
	const core::memory::Scope monitor_group(core::memory::monitors);
	if (!agree(fdtd::app::allocate_planes_and_material_maps(), "cannot allocate the DFT planes / material matrices"))
		return false;
	{
		const core::memory::Scope source_group(core::memory::source);
		if (!agree(fdtd::app::allocate_aux_source(), "cannot allocate the 1D source"))
			return false;
	}

	monitors::custom_boxes().build(_cfg.monitors);

	if (!agree(fdtd::app::allocate_scattering_and_absorption(), "cannot allocate the scattering/absorption monitors"))
		return false;
	if (fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("run", "Memory allocation complete\n");
	return true;
}

void Driver::build_geometry()
{
	auto& c = *_ctx;
	if (fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("run", "Middle layer in y [nm]: %f\n", double(run_cfg.geometry.middle_layer / 1.E-9));

	fdtd::app::create_structures();
	++fdtd::materials::stretch::mesh_generation; // the E scans re-split their rows
	{
		const core::memory::Scope group(core::memory::material_models);
		if (!agree(fdtd::app::allocate_material_state(), "cannot allocate the dispersive material state"))
			return;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// material-index dumps requested by the output settings, the mesh path forces them in DumpGeometry
	if (!mesh_only() && !subdomain.cpml_only)
	{
		write_material_index(c.layout, run_cfg.output.time_plot_3d != 0);
		monitors::custom_boxes().write_masks();
	}
	MPI_Barrier(MPI_COMM_WORLD);

}

void Driver::setup_sources_and_monitors()
{
	auto& c = *_ctx;

	if (run_cfg.materials.magnetic_lorentz)
	{
		fdtd::materials::read_magnetisation_field(run_cfg.output.base_path_input);
		MPI_Barrier(MPI_COMM_WORLD);
		// the static magnetization on the plot planes: materials/magnetization_<plane>.h5
		if (!subdomain.cpml_only)
		{
			if (subdomain.plots_xy())
				save_as_binary_2d(run_cfg.output.materials_dir + "/magnetization_XY.h5", magnetization_x, magnetization_y, magnetization_z,
							   fdtd::grid::output_comm_plot_xy, subdomain.k_animation, c.layout.start, c.layout.global, c.layout.local, c.layout.shift, 0, 1, false, true, true);
			if (subdomain.plots_yz())
				save_as_binary_2d(run_cfg.output.materials_dir + "/magnetization_YZ.h5", magnetization_x, magnetization_y, magnetization_z,
							   fdtd::grid::output_comm_plot_yz, subdomain.i_animation, c.layout.start, c.layout.global, c.layout.local, c.layout.shift, 1, 2, false, true, true);
			if (subdomain.plots_xz())
				save_as_binary_2d(run_cfg.output.materials_dir + "/magnetization_XZ.h5", magnetization_x, magnetization_y, magnetization_z,
							   fdtd::grid::output_comm_plot_xz, subdomain.j_animation, c.layout.start, c.layout.global, c.layout.local, c.layout.shift, 0, 2, false, true, true);
		}
	}

	MPI_Bcast(background_y, run_cfg.grid.nr_cells_uni + 1, MPI_INT, ROOTPROCESS, MPI_COMM_WORLD);

	if (fdtd::grid::world_rank == ROOTPROCESS)
	{
		if (run_cfg.source.tfsf)
			fdtd::io::log("run", "Start FDTD cycle: %s polarization\n", fdtd::app::run_state().cstring);
		if ((run_cfg.source.circ_lh || run_cfg.source.circ_rh) && fdtd::sources::xpol && fdtd::sources::zpol)
			fdtd::io::log("run", "circular polarization: %15.5e s time shift\n", run_cfg.source.circ_time_shift);
		if (run_cfg.source.gauss_beam)
			fdtd::io::log("run", "Start FDTD cycle: Gaussian beam\n");
		else if (run_cfg.source.plane_wave_superposition)
			fdtd::io::log("run", "Start FDTD cycle: Plane wave superposition\n");
		fdtd::io::log("run", "E_peak= %15.5e\n", run_cfg.source.e_peak);

		if (run_cfg.nonlinear.third_harmonic || run_cfg.nonlinear.second_harmonic)
		{
			fdtd::io::log("run", "nonlinear approach (0 = E-field update, 1 = polarization update) = %d\n", run_cfg.nonlinear.pnl_approach);
			fdtd::io::log("run", "meshing (0 = per component, 1 = staircasing) = %d\n", run_cfg.geometry.sta_ir_casing);
			for (int m = 0; m < fdtd::materials::slot_count(); ++m)
			{
				const double* chi = run_cfg.nonlinear.second_harmonic ? fdtd::materials::nonlinear::chi_2[m] : fdtd::materials::nonlinear::chi_3[m];
				fdtd::io::log("run", "material %d %s: nonlinear %d, diagonal only %d, diagonal %15.5e %15.5e %15.5e\n", m, fdtd::materials::slot_name(m).c_str(),
							  int(fdtd::materials::nonlinear::nl[m]), int(run_cfg.nonlinear.second_harmonic ? fdtd::materials::nonlinear::iso_chi_2[m] : fdtd::materials::nonlinear::iso_chi_3[m]),
							  chi[0], chi[7], chi[14]);
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	monitors::time_probes().setup(_cfg.probes);

	fdtd::materials::hydro::nonlocal_barrier_int = 10000;

	if (run_cfg.materials.hydrodynamics || run_cfg.analysis.multipoles)
	{
		if (fdtd::grid::world_rank == ROOTPROCESS)
			fdtd::io::log("run", "Setting up epsilon arrays\n");
		fdtd::materials::hydro::setup_eps_inf_x();
		fdtd::materials::hydro::setup_eps_inf_y();
		fdtd::materials::hydro::setup_eps_inf_z();
		MPI_Barrier(MPI_COMM_WORLD);
		fdtd::grid::exchange_margins(fdtd::materials::hydro::epsinf_x, fdtd::materials::hydro::epsinf_y, fdtd::materials::hydro::epsinf_z);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (run_cfg.thermal.load_temperature)
		fdtd::materials::thermal::load_static_temperature();
	if (run_cfg.thermal.two_temperature_model && run_cfg.thermal.pump_probe && run_cfg.thermal.probe)
		fdtd::materials::thermal::load_probe_temperature_1d();

	if (run_cfg.source.mode_source)
		fdtd::sources::setup_mode_source(); // the MPB profile on the injection plane

	fdtd::materials::hydro::use_backwards_differences = false;
	build_run_lists();
}

// Phase: mesh dump

void Driver::dump_geometry()
{
	auto& c = *_ctx;
	if (!subdomain.cpml_only)
	{
		write_material_index(c.layout, true);
		monitors::custom_boxes().write_masks();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	write_geometry_report();
	fdtd::app::release_run_memory();
}

/// Human-readable qualification report of the mesh: grid, resolution, decomposition, material census.
void Driver::write_geometry_report()
{
	// material census over the interior (TF/SF) cells of every rank, per Yee component
	const int ncodes = fdtd::materials::slot_count();
	std::vector<long long> local(3 * ncodes, 0);
	std::vector<long long> total(3 * ncodes, 0);
	long long out_of_range = 0;
	long long out_of_range_all = 0;
	if (!subdomain.cpml_only)
	{
		int*** comp[3] = {matind_x, matind_y, matind_z};
		for (int a = 0; a < 3; ++a)
			for (int i = 0; i < subdomain.inner_dim_x; ++i)
				for (int j = 0; j < subdomain.inner_dim_y; ++j)
					for (int k = 0; k < subdomain.inner_dim_z; ++k)
					{
						const int m = comp[a][i][j][k];
						if (m >= 0 && m < ncodes)
							++local[a * ncodes + m];
						else
							++out_of_range;
					}
	}
	MPI_Reduce(local.data(), total.data(), 3 * ncodes, MPI_LONG_LONG, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Reduce(&out_of_range, &out_of_range_all, 1, MPI_LONG_LONG, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);

	double rss = peak_rss_mib();
	double rss_max = 0;
	double rss_sum = 0;
	MPI_Reduce(&rss, &rss_max, 1, MPI_DOUBLE, MPI_MAX, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Reduce(&rss, &rss_sum, 1, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);

	// the array tally, per group and in total, last entry
	constexpr double GIB = 1024.0 * 1024.0 * 1024.0;
	const int ng = core::memory::group_count + 1;
	std::vector<unsigned long long> array_bytes(ng, 0);
	std::vector<unsigned long long> array_max(ng, 0);
	std::vector<unsigned long long> array_sum(ng, 0);
	for (int g = 0; g < core::memory::group_count; ++g)
		array_bytes[g] = core::memory::live[g], array_bytes[ng - 1] += core::memory::live[g];
	MPI_Reduce(array_bytes.data(), array_max.data(), ng, MPI_UNSIGNED_LONG_LONG, MPI_MAX, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Reduce(array_bytes.data(), array_sum.data(), ng, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);

	int cpml_only = subdomain.cpml_only ? 1 : 0;
	int cpml_only_all = 0;
	MPI_Reduce(&cpml_only, &cpml_only_all, 1, MPI_INT, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);

	if (fdtd::grid::world_rank != ROOTPROCESS)
		return;

	std::ostringstream r;
	const auto& d = _cfg.grid;
	r << "3D-FDTD geometry report\n\n";
	r << "grid (interior)     : " << d.nx << " x " << d.ny << " x " << d.nz << " cells\n";
	r << "grid (with CPML)    : " << run_cfg.grid.nr_cells_tot_x << " x " << run_cfg.grid.nr_cells_tot_y << " x " << run_cfg.grid.nr_cells_tot_z << " cells\n";
	r << "cell size           : " << d.dx * 1e9 << " x " << d.dy * 1e9 << " x " << d.dz * 1e9 << " nm\n";
	r << "physical size       : " << d.nx * d.dx * 1e9 << " x " << d.ny * d.dy * 1e9 << " x " << d.nz * d.dz * 1e9 << " nm\n";
	r << "time step           : " << run_cfg.time.dt << " s  (courant " << d.courant << ")\n";
	r << "simulated time      : " << run_cfg.time.dt * fdtd::app::run_state().nr_max_time_steps << " s  (" << fdtd::app::run_state().nr_max_time_steps << " steps)\n";

	const double hmax = std::max({d.dx, d.dy, d.dz});
	const double lmin = _cfg.spectrum.lambda_min;
	const double n_host = std::max(1.0, double(fdtd::materials::refractive_index(run_cfg.materials.first_medium, double(PI2 * C0 / lmin)).real()));
	r << "resolution          : " << lmin / hmax << " cells per vacuum wavelength, "
	  << lmin / (n_host * hmax) << " in the first medium (n = " << n_host << ") at lambda_min = " << lmin * 1e9 << " nm\n";
	if (lmin / (n_host * hmax) < 10.0)
		r << "  WARNING: fewer than 10 cells per wavelength in the host medium -> expect numerical dispersion\n";

	r << "\ndecomposition       : " << _cfg.parallel.ranks_x << " x " << _cfg.parallel.ranks_y << " x " << _cfg.parallel.ranks_z << " ranks, " << cpml_only_all << " rank(s) hold CPML only\n";
	{
		char line[128];
		std::snprintf(line, sizeof line, "peak memory (mesh)  : %.3f GiB max per rank, %.3f GiB total\n", rss_max / 1024.0, rss_sum / 1024.0);
		r << line;
	}
	r << "  (the mesh path allocates the same arrays as a simulation, so this is a good estimate of the run;\n"
	     "   sources and monitors built after the mesh are not in it)\n";
	r << "arrays [GiB]        :       largest rank           total\n";
	for (int g = 0; g <= core::memory::group_count; ++g)
	{
		char line[128];
		std::snprintf(line, sizeof line, "  %-18s: %18.3f %15.3f\n", g < core::memory::group_count ? core::memory::GROUP_NAMES[g] : "all arrays",
					  double(array_max[g]) / GIB, double(array_sum[g]) / GIB);
		r << line;
	}
	r << "  in bytes          : " << array_max[ng - 1] << " " << array_sum[ng - 1] << "\n";

	r << "\nmaterial census (interior cells, per Yee component)\n";
	r << "   slot        Ex-cells        Ey-cells        Ez-cells  material, note\n";
	long long structure_cells = 0;
	for (int m = 0; m < ncodes; ++m)
	{
		const long long x = total[m];
		const long long y = total[ncodes + m];
		const long long z = total[2 * ncodes + m];
		if (x + y + z == 0)
			continue;
		std::string note = fdtd::materials::slot_name(m) + (fdtd::materials::nonlinear::nl[m] ? ", nonlinear" : "");
		if (m == run_cfg.materials.antenna_material)
		{
			note += ", structure";
			structure_cells += x + y + z;
		}
		if (m == run_cfg.materials.first_medium)
			note += ", first medium";
		if (m == run_cfg.materials.middle_medium && _cfg.geometry.middle_thickness > 0)
			note += ", middle medium";
		if (m == run_cfg.materials.second_medium)
			note += ", last medium";
		char line[160];
		std::snprintf(line, sizeof(line), "  %5d  %14lld  %14lld  %14lld  %s\n", m, x, y, z, note.c_str());
		r << line;
	}
	if (out_of_range_all)
		r << "  WARNING: " << out_of_range_all << " cell components carry a material slot outside [0, " << ncodes << ")\n";
	if (config::shape_code(_cfg.geometry.shape) != 0 && structure_cells == 0)
		r << "  WARNING: a structure was requested ([geometry] shape = " << _cfg.geometry.shape
		  << ") but no cell carries [geometry] material = " << _cfg.geometry.material
		  << " -> structure not resolved, check dx/dy/dz against the feature sizes\n";
	if (_cfg.multipoles.exact && structure_cells == 0)
		r << "  WARNING: [multipoles] exact = true but the integration volume (structure) is empty\n";

	r << "\nfiles: " << run_cfg.output.materials_dir << "/MaterialIndex*.h5\n";

	const std::string path = run_cfg.output.materials_dir + "/geometry_report.txt";
	std::ofstream(path) << r.str();
	fdtd::io::log("run", "\n%s\nreport written to %s\n", r.str().c_str(), path.c_str());
}

// Phase: finish the run

void Driver::finish_run()
{
	auto& c = *_ctx;
	int j;

	if (run_cfg.thermal.two_temperature_model && run_cfg.thermal.pump_probe && run_cfg.thermal.pump && c.output_pump_data)
	{
		const std::string path = run_cfg.output.base_path_output + "/PumpTemp_" + std::to_string(fdtd::grid::world_rank) + ".txt";
		if (FILE* f = fopen(path.c_str(), "w"))
		{
			for (int ts = 1; ts <= fdtd::app::run_state().nr_max_time_steps; ts++)
				for (j = 0; j <= subdomain.cells_y; j++)
					fprintf(f, "%e\n", fdtd::materials::thermal::temp_save[0][ts][j]);
			fclose(f);
		}
	}

	if (fdtd::app::run_state().ex_file) { fclose(fdtd::app::run_state().ex_file); fdtd::app::run_state().ex_file = NULL; }
	if (fdtd::app::run_state().ez_file) { fclose(fdtd::app::run_state().ez_file); fdtd::app::run_state().ez_file = NULL; }
	MPI_Barrier(MPI_COMM_WORLD);

	// DFT normalisation, then the monitors' analysis and output, collective
	if (!run_cfg.source.dipole_source)
		for (const auto& monitor : c.monitors)
			monitor->normalize();
	MPI_Barrier(MPI_COMM_WORLD);

	fdtd::monitors::count_art = 9; // final ART pass, names the ART-TFSF files
	for (const auto& monitor : c.monitors)
	{
		monitor->finish();
		MPI_Barrier(MPI_COMM_WORLD);
	}

}

// Phase: output + finalize

void Driver::write_outputs()
{
	auto& c = *_ctx;
	for (const auto& monitor : c.monitors)
		monitor->write_final();
	MPI_Barrier(MPI_COMM_WORLD);

	fdtd::app::release_run_memory();
	monitors::time_probes().release();
}

void Driver::finalize()
{
	// a run that reached its last step leaves no state file behind, a later --restart would reuse it
	if (!mesh_only() && !_state_written)
		std::remove((run_cfg.output.state_dir + "/rank_" + std::to_string(_rank) + ".bin").c_str());
	MPI_Barrier(MPI_COMM_WORLD);
	unsigned long long peak = core::memory::peak;
	unsigned long long peak_max = 0;
	unsigned long long peak_sum = 0;
	MPI_Reduce(&peak, &peak_max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Reduce(&peak, &peak_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
	if (fdtd::grid::world_rank == ROOTPROCESS)
	{
		fdtd::io::log("run", "arrays: peak %.3f GiB on the largest rank, %.3f GiB in total (process memory: add the MPI and library overhead per rank)\n",
					  double(peak_max) / 1073741824.0, double(peak_sum) / 1073741824.0);
		const std::string msg = mesh_only() ? std::string("Mesh written to ") + run_cfg.output.materials_dir
						  : _signal ? "Interrupted by signal " + std::to_string(_signal) + " after " + std::to_string(fdtd::app::run_state().timestep) + " iterations; outputs written"
											: "Process " + std::to_string(fdtd::grid::world_rank) + " regularly concluded after " + std::to_string(fdtd::app::run_state().timestep) + " iterations";
		fdtd::io::write_log(msg, run_cfg.output.reports_dir, fdtd::app::run_state().prgstart, "log", true, 1);
		// keep only the directories this run wrote to
		std::error_code ec;
		for (const std::string &dir : {run_cfg.output.reports_dir, run_cfg.output.fields_dir, run_cfg.output.frames_dir, run_cfg.output.materials_dir, run_cfg.output.state_dir})
			if (fs::is_directory(dir, ec) && fs::is_empty(dir, ec))
				fs::remove(dir, ec);
	}
}

} // namespace fdtd
