/**
 * @file app/setup.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @date 2022-06-22
 *
 * @brief Houses functions which allocate the memory for various fields and monitors and the functions that setup and initialize the simulation.
 * All functions allocate the memory to a global address for use externally.
 *  
 *
 * @copyright Copyright 2022
 */
#include "sources/aux_grid.h"
#include "app/run_config.h"
#include "app/run_state.h"
#include "materials/nonlinear/state.h"
#include "io/text.h"
#include "materials/thermal/state.h"
#include "grid/state.h"
#include "materials/state.h"
#include "geometry/state.h"
#include "monitors/state.h"
#include "sources/state.h"
#include "materials/hydro/state.h"
#include "app/setup.h"
#include "sources/signal.h"
#include "materials/materials.h"
#include "materials/spec.h"
#include "sources/spec.h"
#include "grid/lattice.h"
#include "boundary/cpml.h"
#include "monitors/volume_dft.h"
#include "monitors/flux_boxes.h"
#include "monitors/art_planes.h"
#include "fields/fields.h"
#include "grid/subdomain.h"
#include "io/binary.h"
#include "sources/types.h"
#include "app/signal_config.h"
using fdtd::constants::C0;
using fdtd::constants::DEFAULT_TEMPERATURE;
using fdtd::constants::MU0;
using fdtd::constants::PI;
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;
using fdtd::materials::max_slots;
using fdtd::materials::magnetization_x;
using fdtd::materials::magnetization_y;
using fdtd::materials::magnetization_z;
using fdtd::materials::matind;
using fdtd::materials::matind_x;
using fdtd::materials::matind_x_full;
using fdtd::materials::matind_y;
using fdtd::materials::matind_y_full;
using fdtd::materials::matind_z;
using fdtd::materials::matind_z_full;
using fdtd::materials::max_slots;
using fdtd::materials::background_y;

namespace
{
const fdtd::app::RunConfig &run_cfg = fdtd::app::run_config(); ///< the run configuration
} // namespace

namespace
{
fdtd::boundary::Cpml &the_cpml = fdtd::boundary::cpml(); ///< cached: the accessor is out of line
} // namespace

// File-local state, formerly globals used only in this file.
namespace
{
int nr_cells_y_pure;                                        ///< Global cells along y without the CPML.
std::vector<int> deferred_poles;                            ///< materials whose poles allocate_material_state allocates
int blocklength;                                            ///< Scratch: block length for the MPI_Type_vector calls in setup.
int stride;                                                 ///< Scratch: stride for the MPI_Type_vector calls in setup.
} // namespace

 namespace
{
/// state of dispersive slot mm stored for the rows of its cells only: the point-update models,
/// Drude + critical points, Drude, Lorentz, modified Lorentz, without the nonlocal / hydrodynamic /
/// two-temperature arrays, which take neighbouring cells
bool boxed_state(int mm)
{
	const auto &rc = fdtd::app::run_config();
	return fdtd::materials::kind(mm) < fdtd::materials::Kind::hydrodynamic && !rc.materials.nonlocal && !rc.thermal.two_temperature_model && !rc.materials.hydrodynamics;
}

/// the slot carries a free-electron current J: hydrodynamic, Drude-Lorentz
bool electron_current(int m)
{
	return fdtd::materials::is(m, fdtd::materials::Kind::hydrodynamic) || fdtd::materials::is(m, fdtd::materials::Kind::drude_lorentz);
}

/// the slot carries pole polarisations: every dispersive model but Drude and Lorentz
bool pole_state(int m)
{
	return !fdtd::materials::is(m, fdtd::materials::Kind::drude) && !fdtd::materials::is(m, fdtd::materials::Kind::lorentz);
}
} // namespace

/**
  * @brief It's executed by every process to initialize all local parameters for the simulation.
  *
  * @return true when initialization is complete
  * @return false when error occurs during intialization
  * 
  * @details Completes the run configuration, derived flags, steps, cell sizes, and sets this rank's
  * local values.
  */
bool fdtd::app::initialise_run()
{
	fdtd::app::RunConfig &run_cfg = fdtd::app::edit_run_config(); // configure phase
	// INIZIALIZATION - START

	subdomain.cells_x = fdtd::grid::cells_per_rank_x[subdomain.position[0]];
	subdomain.cells_y = fdtd::grid::cells_per_rank_y[subdomain.position[1]];
	subdomain.cells_z = fdtd::grid::cells_per_rank_z[subdomain.position[2]];

	// Flags for source injection types
	run_cfg.source.tfsf = run_cfg.source.kind == fdtd::sources::kind::plane_wave;
	run_cfg.source.dipole_source = run_cfg.source.kind == fdtd::sources::kind::dipole;
	run_cfg.source.gauss_beam = run_cfg.source.kind == fdtd::sources::kind::gaussian_beam;
	run_cfg.source.plane_wave_superposition = run_cfg.source.kind == fdtd::sources::kind::wave_superposition;

	// activate 1D simulations if settings are given, or if enforced using a key in pphinfoini.json
	run_cfg.source.oblique_incidence_with_medium = (run_cfg.source.plane_wave_superposition && ( run_cfg.materials.middle_medium != run_cfg.materials.first_medium || run_cfg.materials.second_medium != run_cfg.materials.first_medium)) || run_cfg.source.enforce_one_dim_simulation;

	run_cfg.analysis.dft_plot = run_cfg.analysis.dft_plane_mask != 0;
	run_cfg.nonlinear.second_harmonic = run_cfg.nonlinear.shg != 0;
	run_cfg.nonlinear.third_harmonic = run_cfg.nonlinear.thg != 0;

	{
		fdtd::materials::LoadSpec spec;
		spec.harmonic_order = run_cfg.nonlinear.second_harmonic ? 2 : run_cfg.nonlinear.third_harmonic ? 3 : 0;
		spec.dy = fdtd::grid::lattice().d[1];
		spec.fermi_velocity_scale = run_cfg.materials.fermi_v;
		fdtd::materials::load_slots(run_cfg.materials.definitions, spec);
	}

	run_cfg.materials.magnetic_lorentz = fdtd::materials::is(run_cfg.materials.antenna_material, fdtd::materials::Kind::magnetic_lorentz);

	if (run_cfg.materials.magnetic_lorentz && fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("setup", "Magnetic Lorentz model\n");

	run_cfg.geometry.grating = false;					// Legacy

	if (run_cfg.geometry.antenna_type == fdtd::geometry::shape::round_hybrid && run_cfg.materials.antenna_material == run_cfg.materials.second_medium && run_cfg.grid.pbc_x && run_cfg.grid.pbc_z)
		run_cfg.geometry.grating = true;

	fdtd::monitors::single_photon = false; // was the bullseye shape, removed

	MPI_Bcast(&run_cfg.source.oblique_incidence_with_medium, 1, MPI_C_BOOL, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&run_cfg.source.center_freq, 1, MPI_DOUBLE, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&run_cfg.source.pulse_width, 1, MPI_DOUBLE, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&run_cfg.source.chirp, 1, MPI_DOUBLE, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&run_cfg.source.chirp_2, 1, MPI_DOUBLE, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&fdtd::monitors::sphere_rad_multipole, 1, MPI_INT, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&fdtd::monitors::nphi_multipole, 1, MPI_INT, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&fdtd::monitors::ntheta_multipole, 1, MPI_INT, ROOTPROCESS, MPI_COMM_WORLD);
	//exact multipoles
	MPI_Bcast(&fdtd::monitors::exact_multipoles, 1, MPI_C_BOOL, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&fdtd::monitors::lwa_multipoles, 1, MPI_C_BOOL, ROOTPROCESS, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	fdtd::monitors::phase_normalization = run_cfg.output.dft_phase_norm != 0;

	if (run_cfg.grid.pbc_x || run_cfg.grid.pbc_z)
		run_cfg.analysis.art = true; //only on the plane to calculate R and T
	else
		run_cfg.analysis.art = false;
	if (!run_cfg.grid.pbc_x && !run_cfg.grid.pbc_z && run_cfg.source.kind != fdtd::sources::kind::mode && (run_cfg.source.kind != fdtd::sources::kind::none || (run_cfg.source.kind == fdtd::sources::kind::none && run_cfg.geometry.antenna_type == fdtd::geometry::shape::round_hybrid)))
		fdtd::monitors::q_coeff = true; // se true calcola DFT ovunque --> da modificare!&& !SHG && !THG
	else
		fdtd::monitors::q_coeff = false;
	fdtd::monitors::dft_plane_xz = false; //only one plane at the end of the FDTD iterations
	fdtd::monitors::dft_plane_yz = false;
	fdtd::monitors::dft_plane_xy = false;
	fdtd::monitors::dft_xz_1 = false;
	fdtd::monitors::dft_xz_2 = false;
	fdtd::monitors::dft_h = false;
	if (run_cfg.analysis.dft_plot)
	{
		// bitmask from [output] dft_yz / dft_xz / dft_xy / dft_magnetic, see DftPlaneMask
		fdtd::monitors::dft_plane_yz = (run_cfg.analysis.dft_plane_mask & 1) != 0;
		fdtd::monitors::dft_plane_xz = (run_cfg.analysis.dft_plane_mask & 2) != 0;
		fdtd::monitors::dft_plane_xy = (run_cfg.analysis.dft_plane_mask & 4) != 0;
		fdtd::monitors::dft_h = (run_cfg.analysis.dft_plane_mask & 8) != 0;
	}

	if (fdtd::monitors::dft_plane_xz && run_cfg.geometry.is(fdtd::geometry::shape::cylinder))
	{
		fdtd::monitors::dft_xz_1 = true;
		fdtd::monitors::dft_xz_2 = true;
	}
	if (fdtd::monitors::dft_plane_xz && (run_cfg.materials.first_medium != run_cfg.materials.middle_medium && run_cfg.materials.middle_medium != run_cfg.materials.second_medium) && run_cfg.geometry.is(fdtd::geometry::shape::round_hybrid))
	{
		fdtd::monitors::dft_xz_1 = true;
	}

	if (fdtd::monitors::dft_plane_xz && (run_cfg.geometry.is(fdtd::geometry::shape::rounded_monopole) || (run_cfg.geometry.is(fdtd::geometry::shape::rounded_dipole) || run_cfg.geometry.is(fdtd::geometry::shape::rounded_dipole_offset_planes)) || run_cfg.geometry.is(fdtd::geometry::shape::rounded_triangle)))
	{
		fdtd::monitors::dft_xz_1 = true;
		fdtd::monitors::dft_xz_2 = true;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (run_cfg.source.tfsf)
	{
		if (run_cfg.source.polarization_deg == double(0))
		{
			fdtd::sources::xpol = true;
			fdtd::sources::zpol = false;
		}
		if (run_cfg.source.polarization_deg == double(90))
		{
			fdtd::sources::xpol = false;
			fdtd::sources::zpol = true;
		}
		if (run_cfg.source.polarization_deg != double(90) && run_cfg.source.polarization_deg != double(0))
		{
			fdtd::sources::xpol = true;
			fdtd::sources::zpol = true;
			fdtd::sources::polarization_angle = double(run_cfg.source.polarization_deg * PI / 180);
		}
	}
	fdtd::materials::thermal::str_outfilepathname_t = run_cfg.output.state_dir + "/electron_temperature_";
	fdtd::materials::thermal::str_infilepathname_t = run_cfg.output.base_path_input + "/electron_temperature";
	fdtd::materials::thermal::str_fileextension = ".h5";

	fdtd::materials::laplacian = false; // This should always be false, but I'd like to keep the booleans just in case
	fdtd::materials::gradient_divergence = true; // This should always be true, but I'd like to keep the booleans just in case


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	set_time_step();
	{
		fdtd::materials::ModelSpec ms;
		ms.dt = run_cfg.time.dt;
		ms.fermi_velocity = run_cfg.materials.fermi_v;
		ms.shg = run_cfg.nonlinear.second_harmonic, ms.thg = run_cfg.nonlinear.third_harmonic;
		ms.harmonic = run_cfg.nonlinear.seperate_linear_nonlinear;
		ms.cellwise = run_cfg.source.mode_source;
		ms.nonlocal = run_cfg.materials.nonlocal, ms.hydrodynamics = run_cfg.materials.hydrodynamics, ms.gnor = run_cfg.materials.gnor;
		ms.two_temperature = run_cfg.thermal.two_temperature_model, ms.load_temperature = run_cfg.thermal.load_temperature;
		ms.pump = run_cfg.thermal.pump, ms.probe = run_cfg.thermal.probe, ms.pump_probe = run_cfg.thermal.pump_probe;
		ms.antenna_material = run_cfg.materials.antenna_material;
		ms.medium_first = run_cfg.materials.first_medium, ms.medium_second = run_cfg.materials.second_medium;
		ms.staircase = run_cfg.geometry.sta_ir_casing;
		ms.tfsf_y_shrink = run_cfg.source.tfsf_y_shrink;
		fdtd::materials::set_model_spec(ms);
	}
	fdtd::app::configure_signal(); // after the time step: the signal keeps dt

	if (run_cfg.thermal.two_temperature_model && run_cfg.thermal.pump && run_cfg.thermal.pump_probe)
	{
		fdtd::app::run_state().nr_max_time_steps += int(220e-15 / run_cfg.time.dt);
		fdtd::io::log("setup", "Increasing Time Steps to %d\n", fdtd::app::run_state().nr_max_time_steps);
	}
	run_cfg.time.coeff_h = run_cfg.time.dt / MU0;
	MPI_Barrier(MPI_COMM_WORLD);
	// origin, in number of cells, of the local process with respect to the global reference system
	subdomain.start_x = 0;
	subdomain.start_y = 0;
	subdomain.start_z = 0;
	for (int m = 0; m < subdomain.position[0]; ++m)
		subdomain.start_x += fdtd::grid::cells_per_rank_x[m];
	for (int m = 0; m < subdomain.position[1]; ++m)
		subdomain.start_y += fdtd::grid::cells_per_rank_y[m];
	for (int m = 0; m < subdomain.position[2]; ++m)
		subdomain.start_z += fdtd::grid::cells_per_rank_z[m];
	// physical origin of the local process with respect to the global reference system
	subdomain.origin_x = subdomain.start_x * run_cfg.grid.dx;
	subdomain.origin_y = subdomain.start_y * run_cfg.grid.dy;
	subdomain.origin_z = subdomain.start_z * run_cfg.grid.dz;

	// physical size of the local process
	subdomain.length_x = subdomain.cells_x * run_cfg.grid.dx;
	subdomain.length_y = subdomain.cells_y * run_cfg.grid.dy;
	subdomain.length_z = subdomain.cells_z * run_cfg.grid.dz;
	subdomain.end_x = subdomain.start_x + subdomain.cells_x;
	subdomain.end_y = subdomain.start_y + subdomain.cells_y;
	subdomain.end_z = subdomain.start_z + subdomain.cells_z;

	// communicators for parallel output: the ranks that own the plane/face, decided by each rank for itself
	fdtd::grid::output_comm = fdtd::grid::member_comm(subdomain.writes_output());
	fdtd::grid::output_comm_plot_yz = fdtd::grid::member_comm(subdomain.plots_yz());
	fdtd::grid::output_comm_plot_xz = fdtd::grid::member_comm(subdomain.plots_xz());
	fdtd::monitors::output_comm_plot_xz_1 = fdtd::grid::member_comm(subdomain.plots_xz_1());
	fdtd::grid::output_comm_plot_xz_2 = fdtd::grid::member_comm(subdomain.plots_xz_2());
	fdtd::grid::output_comm_plot_xy = fdtd::grid::member_comm(subdomain.plots_xy());
	fdtd::monitors::output_comm_sca_nea_x = fdtd::grid::member_comm(subdomain.flux_sca_nea_x());
	fdtd::monitors::output_comm_sca_far_x = fdtd::grid::member_comm(subdomain.flux_sca_far_x());
	fdtd::monitors::output_comm_sca_nea_y = fdtd::grid::member_comm(subdomain.flux_sca_nea_y());
	fdtd::monitors::output_comm_sca_far_y = fdtd::grid::member_comm(subdomain.flux_sca_far_y());
	fdtd::monitors::output_comm_sca_nea_z = fdtd::grid::member_comm(subdomain.flux_sca_nea_z());
	fdtd::monitors::output_comm_sca_far_z = fdtd::grid::member_comm(subdomain.flux_sca_far_z());
	if (!subdomain.cpml_only)
	{
		subdomain.output_size[0] = subdomain.inner_dim_x + 1;
		subdomain.output_size[1] = subdomain.inner_dim_y + 1;
		subdomain.output_size[2] = subdomain.inner_dim_z + 1;
		subdomain.output_start[0] = subdomain.shift_output_x;
		subdomain.output_start[1] = subdomain.shift_output_y;
		subdomain.output_start[2] = subdomain.shift_output_z;

	}
	// IT DEFINES THE TYPES FOR THE EXCHANGE OF DATA BETWEEN PROCESSES
	int vector_count; // elements of each halo-plane MPI_Type_vector below
	vector_count = (subdomain.cells_x + 1) * (subdomain.cells_y + 1);
	blocklength = 1;
	stride = subdomain.cells_z + 1;
	MPI_Type_vector(vector_count, blocklength, stride, MPI_DOUBLE, &fdtd::grid::xy_interface_data_type);
	MPI_Type_commit(&fdtd::grid::xy_interface_data_type);

	vector_count = (subdomain.cells_x + 3) * (subdomain.cells_y + 3);
	blocklength = 1;
	stride = subdomain.cells_z + 3;
	MPI_Type_vector(vector_count, blocklength, stride, MPI_DOUBLE, &fdtd::grid::xy_interface_data_type_2);
	MPI_Type_commit(&fdtd::grid::xy_interface_data_type_2);

	vector_count = subdomain.cells_x + 1;
	blocklength = subdomain.cells_z + 1;
	stride = (subdomain.cells_y + 1) * (subdomain.cells_z + 1);
	MPI_Type_vector(vector_count, blocklength, stride, MPI_DOUBLE, &fdtd::grid::xz_interface_data_type);
	MPI_Type_commit(&fdtd::grid::xz_interface_data_type);

	vector_count = subdomain.cells_x + 3;
	blocklength = subdomain.cells_z + 3;
	stride = (subdomain.cells_y + 3) * (subdomain.cells_z + 3);
	MPI_Type_vector(vector_count, blocklength, stride, MPI_DOUBLE, &fdtd::grid::xz_interface_data_type_2);
	MPI_Type_commit(&fdtd::grid::xz_interface_data_type_2);

	vector_count = (subdomain.cells_y + 1) * (subdomain.cells_z + 1);
	blocklength = 1;
	stride = 1;
	MPI_Type_vector(vector_count, blocklength, stride, MPI_DOUBLE, &fdtd::grid::yz_interface_data_type);
	MPI_Type_commit(&fdtd::grid::yz_interface_data_type);

	vector_count = (subdomain.cells_y + 3) * (subdomain.cells_z + 3);
	blocklength = 1;
	stride = 1;
	MPI_Type_vector(vector_count, blocklength, stride, MPI_DOUBLE, &fdtd::grid::yz_interface_data_type_2);
	MPI_Type_commit(&fdtd::grid::yz_interface_data_type_2);

	fdtd::grid::neighbour_ranks(subdomain.position, subdomain.neighbours_x, subdomain.neighbours_y, subdomain.neighbours_z);
	fdtd::grid::neighbour_ranks_nl(subdomain.position, subdomain.neighbours_nl_x, subdomain.neighbours_nl_y, subdomain.neighbours_nl_z);
	// [02]

	fdtd::materials::compute_ade_coefficients();
	{
		const double dt = run_cfg.time.dt; // the TF/SF H corrections, were set with the material coefficients
		const fdtd::grid::Lattice &lat = fdtd::grid::lattice();
		fdtd::sources::chx = dt / MU0 / lat.d[0];
		fdtd::sources::chy = dt / MU0 / lat.d[1];
		fdtd::sources::chz = dt / MU0 / lat.d[2];
	}

	if (fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("setup", "dimensions without CPML: %d %d %d\n", fdtd::grid::lattice().inner(0), fdtd::grid::lattice().inner(1), fdtd::grid::lattice().inner(2));
	if (fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("setup", "dimensions CPML: %d %d %d\n", run_cfg.grid.num_cpml_cells_x, run_cfg.grid.num_cpml_cells_y, run_cfg.grid.num_cpml_cells_z);
	if (fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("setup", "dimensions PW: %d %d %d\n", run_cfg.grid.nr_inflate_x, run_cfg.grid.nr_inflate_y, run_cfg.grid.nr_inflate_z);
	MPI_Barrier(MPI_COMM_WORLD);

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	subdomain.owned_cells_x = subdomain.cells_x;
	if (subdomain.position[0] == fdtd::grid::num_mpi_processes_x - 1)
		subdomain.owned_cells_x = subdomain.cells_x - 1;

	subdomain.owned_cells_y = subdomain.cells_y;
	if (subdomain.position[1] == fdtd::grid::num_mpi_processes_y - 1)
		subdomain.owned_cells_y = subdomain.cells_y - 1;

	subdomain.owned_cells_z = subdomain.cells_z;
	if (subdomain.position[2] == fdtd::grid::num_mpi_processes_z - 1)
		subdomain.owned_cells_z = subdomain.cells_z - 1;

	subdomain.has_plane_wave = false;
	nr_cells_y_pure = run_cfg.grid.nr_cells_tot_y - 2 * run_cfg.grid.num_cpml_cells_y;
	run_cfg.grid.nr_cells_uni = run_cfg.source.pw_factor * run_cfg.grid.nr_cells_tot_y; // pw_factor*nr_cells_y_pure;
	fdtd::sources::nr_cell_uni = run_cfg.grid.nr_cells_uni - 1;
	fdtd::sources::nr_cell_par_uni = run_cfg.grid.nr_cells_uni - run_cfg.grid.num_cpml_cells_y;
	fdtd::sources::j_min_cpml_uni = run_cfg.grid.num_cpml_cells_y;
	fdtd::sources::j_max_cpml_uni = run_cfg.grid.nr_cells_uni - run_cfg.grid.num_cpml_cells_y;
	fdtd::geometry::j_offset_uni = int(run_cfg.grid.nr_cells_uni / 2) - int(nr_cells_y_pure / 2);
	fdtd::sources::jpw_start_uni = fdtd::geometry::j_offset_uni + run_cfg.grid.nr_inflate_y;
	fdtd::sources::jpw_source_uni = fdtd::sources::jpw_start_uni - 1;
	fdtd::sources::jpw_end_uni = fdtd::geometry::j_offset_uni + nr_cells_y_pure - run_cfg.grid.nr_inflate_y - run_cfg.source.tfsf_y_shrink;
	subdomain.j_for_uni = subdomain.jpw_pad - run_cfg.grid.num_cpml_cells_y + fdtd::geometry::j_offset_uni;

	if ((subdomain.pw_nea_x || subdomain.pw_far_x || subdomain.pw_nea_y || subdomain.pw_far_y || subdomain.pw_nea_z || subdomain.pw_far_z) && !run_cfg.source.dipole_source)
	{
		subdomain.has_plane_wave = true;

		fdtd::sources::ipw_near_h = subdomain.ipw_1 - 1;
		fdtd::sources::jpw_near_h = subdomain.jpw_1 - 1;
		fdtd::sources::kpw_near_h = subdomain.kpw_1 - 1;
		fdtd::sources::ipw_min_mod = subdomain.ipw_min;
		if (subdomain.ipw_min == 0)
			fdtd::sources::ipw_min_mod += 1;
		fdtd::sources::kpw_min_mod = subdomain.kpw_min;
		if (subdomain.kpw_min == 0)
			fdtd::sources::kpw_min_mod += 1;
	}

	if (subdomain.qsca_nea_x || subdomain.qsca_far_x || subdomain.qsca_nea_y || subdomain.qsca_far_y || subdomain.qsca_nea_z || subdomain.qsca_far_z)
	{
		//for near-to-far surface
		fdtd::monitors::globalsizes_nf[0] = run_cfg.grid.nr_cells_tot_x - 2 * run_cfg.grid.num_cpml_cells_x - 2 * run_cfg.analysis.nr_qsca_x + 1;
		fdtd::monitors::globalsizes_nf[1] = run_cfg.grid.nr_cells_tot_y - 2 * run_cfg.grid.num_cpml_cells_y - 2 * run_cfg.analysis.nr_qsca_y + 1;
		fdtd::monitors::globalsizes_nf[2] = run_cfg.grid.nr_cells_tot_z - 2 * run_cfg.grid.num_cpml_cells_z - 2 * run_cfg.analysis.nr_qsca_z + 1;
		fdtd::monitors::localsize_nf[0] = subdomain.i_qsca_max - subdomain.i_qsca_min + 1; //subdomain.inner_dim_x+1;
		fdtd::monitors::localsize_nf[1] = subdomain.j_qsca_max - subdomain.j_qsca_min + 1; //subdomain.inner_dim_y+1;
		fdtd::monitors::localsize_nf[2] = subdomain.k_qsca_max - subdomain.k_qsca_min + 1; //subdomain.inner_dim_z+1;
		fdtd::monitors::startindices_nf[0] = subdomain.shift_output_x_nf;
		fdtd::monitors::startindices_nf[1] = subdomain.shift_output_y_nf;
		fdtd::monitors::startindices_nf[2] = subdomain.shift_output_z_nf;

		//if(subdomain.qsca_nea_x && (  (subdomain.j_qsca_min != subdomain.j_qsca_max) &&
		//  subdomain.k_qsca_min != subdomain.k_qsca_max))
	}

	{
		fdtd::sources::Spec sp;
		const auto &src = run_cfg.source;
		sp.dt = run_cfg.time.dt, sp.coeff_h = run_cfg.time.coeff_h;
		sp.one_dim_steps = run_cfg.time.one_dim_nr_max_time_steps, sp.one_dim_cells = run_cfg.grid.nr_cells_uni;
		sp.medium[0] = run_cfg.materials.first_medium, sp.medium[1] = run_cfg.materials.middle_medium, sp.medium[2] = run_cfg.materials.second_medium;
		sp.middle_layer = run_cfg.geometry.middle_layer;
		sp.centre[0] = run_cfg.geometry.center_position_x, sp.centre[1] = run_cfg.geometry.center_position_y, sp.centre[2] = run_cfg.geometry.center_position_z;
		sp.grating = run_cfg.geometry.grating;
		sp.hydrodynamics = run_cfg.materials.hydrodynamics;
		sp.harmonic = run_cfg.nonlinear.seperate_linear_nonlinear;
		sp.shg = run_cfg.nonlinear.second_harmonic, sp.thg = run_cfg.nonlinear.third_harmonic;
		sp.nr_frequencies = run_cfg.analysis.nr_an_frequencies;
		sp.fundamental_frequency = src.fundamental_frequency, sp.highest_frequency = src.highest_frequency, sp.start_cw = src.start_cw;
		sp.tfsf = src.tfsf, sp.mode_source = src.mode_source, sp.dipole = src.dipole_source, sp.gauss_beam = src.gauss_beam;
		sp.superposition = src.plane_wave_superposition, sp.oblique_with_medium = src.oblique_incidence_with_medium;
		sp.circ_time_shift = src.circ_time_shift, sp.signal_type = src.signal_type;
		sp.cw_pulse = src.cw_pulse, sp.continuous_wave = src.continuous_wave;
		sp.theta_intervals = src.theta_nr_intervals, sp.phi_intervals = src.phi_nr_intervals;
		sp.theta_max_degrees = src.theta_max_degrees, sp.phi_max_degrees = src.phi_max_degrees;
		sp.waist_factor = src.beam_w_0_coeff, sp.focal_factor = src.focal_length_coeff;
		sp.plot_xy = run_cfg.output.plot_xy, sp.plot_xz = run_cfg.output.plot_xz, sp.plot_yz = run_cfg.output.plot_yz;
		sp.reports_dir = run_cfg.output.reports_dir, sp.input_dir = run_cfg.output.base_path_input;
		sp.log_start = fdtd::app::run_state().prgstart;
		fdtd::sources::set_spec(sp);
	}
	return true;
}
/**
 * @brief It sets the simulation up based on the structure to simulate, valid only for some structures
 * @details Based on the type of structure (run_cfg.geometry.antenna_type), for example photonic crystals,
 * some simulation parameters, unit cell size, type of excitation, can be set up automatically.
 */
void fdtd::app::setup_domain()
{
	fdtd::app::RunConfig &run_cfg = fdtd::app::edit_run_config(); // configure phase
	run_cfg.source.tfsf_y_shrink = 0;
	if (run_cfg.source.kind > 10)
	{
		run_cfg.source.tfsf_y_shrink = run_cfg.source.kind;
		run_cfg.source.kind = 1;
	}

//automatic CPML setup

	if (run_cfg.source.kind == fdtd::sources::kind::plane_wave && run_cfg.grid.pbc_x == 1)
	{
		run_cfg.grid.nr_cells_tot_x = run_cfg.grid.nr_cells_tot_x - 2 * run_cfg.grid.num_cpml_cells_x + 2;
		run_cfg.grid.num_cpml_cells_x = 1;
		run_cfg.grid.nr_inflate_x = 0;
	}
	if (run_cfg.source.kind == fdtd::sources::kind::plane_wave && run_cfg.grid.pbc_z == 1)
	{
		run_cfg.grid.nr_cells_tot_z = run_cfg.grid.nr_cells_tot_z - 2 * run_cfg.grid.num_cpml_cells_z + 2;
		run_cfg.grid.num_cpml_cells_z = 1;
		run_cfg.grid.nr_inflate_z = 0;
	}
}

 /**
  * @brief Allocates the memory for the fields for different kind of models and simulations
  *
  * @return true when memory allocation has been succesfull
  * @return false when bad_alloc is thrown
  */
bool fdtd::app::allocate_fields()
{
	auto &ds = fdtd::materials::dispersive_state();
	const auto &mp = fdtd::materials::model_parameters();
	int m;
	int n;
	int i;
	int j;
	int k;
	int checker = 0;

	try
	{
		fdtd::fields::Needs needs; // E/H, their history and scratch copies: fields module
		needs.harmonic = run_cfg.nonlinear.seperate_linear_nonlinear;
		needs.e_prev_2 = run_cfg.materials.magnetic_lorentz;
		needs.h_prev = run_cfg.materials.hydrodynamics;
		needs.e_scratch = run_cfg.nonlinear.second_harmonic || run_cfg.nonlinear.third_harmonic || run_cfg.materials.magnetic_lorentz;
		{
			const fdtd::core::memory::Scope group(fdtd::core::memory::fields);
			if (!fdtd::fields::fields().allocate(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z, needs))
				return false;
		}
		const fdtd::core::memory::Scope group(fdtd::core::memory::material_models); // the rest of this function, but the maps

		//-Izzat-

		if (run_cfg.materials.magnetic_lorentz)
		{
			magnetization_x = fdtd::field3<double>(subdomain.cells_x + 1, subdomain.cells_y + 1, subdomain.cells_z + 1);
			magnetization_y = fdtd::field3<double>(subdomain.cells_x + 1, subdomain.cells_y + 1, subdomain.cells_z + 1);
			magnetization_z = fdtd::field3<double>(subdomain.cells_x + 1, subdomain.cells_y + 1, subdomain.cells_z + 1);
		}

		{
			const fdtd::core::memory::Scope maps(fdtd::core::memory::material_maps);
			fdtd::materials::material_map().allocate_full(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
		}

		if (run_cfg.thermal.two_temperature_model || run_cfg.thermal.load_temperature)
		{

			if (run_cfg.thermal.pump_probe)
			{
				if (run_cfg.thermal.pump)
					fdtd::materials::thermal::temp_save = fdtd::field3<double>(1, fdtd::app::run_state().nr_max_time_steps, subdomain.cells_y);
				fdtd::materials::thermal::nr_max_time_steps_pump = fdtd::app::run_state().nr_max_time_steps + int(220e-15 / run_cfg.time.dt);
				if (run_cfg.thermal.probe)
					fdtd::materials::thermal::temp_save = fdtd::field3<double>(1, fdtd::materials::thermal::nr_max_time_steps_pump, subdomain.cells_y);
			}

			fdtd::materials::thermal::electron_temp = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
			fdtd::materials::thermal::lattice_temp = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
			fdtd::monitors::energy_density = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
			if (!run_cfg.materials.hydrodynamics)
			{
				fdtd::materials::hydro::ex_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
				fdtd::materials::hydro::ey_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
				fdtd::materials::hydro::ez_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			}
			for (i = 0; i <= subdomain.cells_x; i++)
			{
				for (j = 0; j <= subdomain.cells_y; j++)
				{
					for (k = 0; k <= subdomain.cells_z; k++)
					{
						fdtd::materials::thermal::electron_temp[i][j][k] = DEFAULT_TEMPERATURE;
						fdtd::materials::thermal::lattice_temp[i][j][k] = DEFAULT_TEMPERATURE;
					}
				}
			}
		}

		if (!run_cfg.materials.hydrodynamics && run_cfg.analysis.multipoles && !run_cfg.thermal.two_temperature_model)
		{
			fdtd::materials::hydro::ex_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::ey_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::ez_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);

			fdtd::materials::hydro::px_cp_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::py_cp_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::pz_cp_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);

			fdtd::materials::hydro::epsinf_x = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::epsinf_y = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::epsinf_z = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);

			fdtd::materials::hydro::nd_x = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
			fdtd::materials::hydro::nd_y = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
			fdtd::materials::hydro::nd_z = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
		}

		if (run_cfg.materials.hydrodynamics)
		{

			fdtd::materials::hydro::ex_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::ey_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::ez_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);

			fdtd::materials::hydro::px_cp_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::py_cp_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::pz_cp_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);

			if (run_cfg.nonlinear.seperate_linear_nonlinear)
			{
				fdtd::materials::hydro::px_cp_2_nl = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
				fdtd::materials::hydro::py_cp_2_nl = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
				fdtd::materials::hydro::pz_cp_2_nl = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);

				fdtd::materials::hydro::ex_2_nl = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
				fdtd::materials::hydro::ey_2_nl = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
				fdtd::materials::hydro::ez_2_nl = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			}

			fdtd::materials::hydro::epsinf_x = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::epsinf_y = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
			fdtd::materials::hydro::epsinf_z = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);
		}

		for (m = 0; m < fdtd::materials::slot_count(); ++m)
		{
			if (fdtd::materials::dispersive(m))
			{
				if (!run_cfg.materials.nonlocal && !run_cfg.thermal.two_temperature_model && !run_cfg.materials.hydrodynamics)
				{
					ds.c[0].p.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
					ds.c[1].p.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
					ds.c[2].p.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
					if (electron_current(m))
					{

						ds.c[0].j.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
						ds.c[1].j.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
						ds.c[2].j.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//

						fdtd::materials::hydro::jx_prev.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
						fdtd::materials::hydro::jy_prev.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
						fdtd::materials::hydro::jz_prev.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
					}
					if (fdtd::materials::is(m, fdtd::materials::Kind::magnetic_lorentz))
					{
						ds.c[0].p_prev_2.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
						ds.c[1].p_prev_2.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
						ds.c[2].p_prev_2.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);

						ds.c[0].p_prev_3.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
						ds.c[1].p_prev_3.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
						ds.c[2].p_prev_3.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
					}
				}
				else
				{ //if nonlocal
					if (electron_current(m))
					{
						if (run_cfg.materials.hydrodynamics || run_cfg.thermal.two_temperature_model)
						{
							ds.c[0].j.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
							ds.c[2].j.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
							ds.c[1].j.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//

							fdtd::materials::hydro::jx_prev.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
							fdtd::materials::hydro::jy_prev.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
							fdtd::materials::hydro::jz_prev.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//

							fdtd::materials::hydro::jx_extra.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
							fdtd::materials::hydro::jy_extra.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
							fdtd::materials::hydro::jz_extra.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//

							if (run_cfg.nonlinear.seperate_linear_nonlinear)
							{
								ds.c[0].j_nl.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
								ds.c[2].j_nl.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
								ds.c[1].j_nl.require(m, subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//

							}
						}

						else
						{
							ds.c[0].j.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
							ds.c[2].j.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
							ds.c[1].j.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//

							fdtd::materials::hydro::jx_prev.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
							fdtd::materials::hydro::jy_prev.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
							fdtd::materials::hydro::jz_prev.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
						}
					}
					else
					{
						ds.c[0].p.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
						ds.c[1].p.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); //subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
						ds.c[2].p.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
					}

					if (checker == 0)
					{
						if (!run_cfg.materials.hydrodynamics && !run_cfg.thermal.two_temperature_model)
						{
							fdtd::materials::hydro::px_nl = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
							fdtd::materials::hydro::py_nl = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); //subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
							fdtd::materials::hydro::pz_nl = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
							ds.c[0].p_nl_prev = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
							ds.c[1].p_nl_prev = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); //subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
							ds.c[2].p_nl_prev = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
							fdtd::materials::hydro::px_nl_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
							fdtd::materials::hydro::py_nl_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); //subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
							fdtd::materials::hydro::pz_nl_2 = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
						}

						if (run_cfg.materials.hydrodynamics)
						{																				  //Electron number density
							fdtd::materials::hydro::nd_x = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
							fdtd::materials::hydro::nd_y = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
							fdtd::materials::hydro::nd_z = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2);	   // subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//
							fdtd::materials::hydro::nd_z_prev = fdtd::field3<double>(subdomain.cells_x + 2, subdomain.cells_y + 2, subdomain.cells_z + 2); //subdomain.inner_dim_x+2, subdomain.inner_dim_y+2, subdomain.inner_dim_z+2);//

							//REMOVE

						}

						if (run_cfg.materials.gnor)
						{
							fdtd::materials::hydro::grad_div_x = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
							fdtd::materials::hydro::grad_div_y = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
							fdtd::materials::hydro::grad_div_z = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);

							fdtd::materials::hydro::grad_div_x_n = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
							fdtd::materials::hydro::grad_div_y_n = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
							fdtd::materials::hydro::grad_div_z_n = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);

							fdtd::materials::hydro::grad_div_x_n_1 = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
							fdtd::materials::hydro::grad_div_y_n_1 = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
							fdtd::materials::hydro::grad_div_z_n_1 = fdtd::field3<double>(subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);
						}
						checker = 1;
					}
				}
				ds.c[0].p_prev.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
				ds.c[1].p_prev.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//
				ds.c[2].p_prev.require(m, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);//

				if (pole_state(m) && boxed_state(m))
					deferred_poles.push_back(m); // allocated after the geometry, rows of the material only
				else if (pole_state(m))
				{
					for (n = 0; n < mp.nr_cp[m]; ++n)
					{
						if (!ds.c[0].pc[m]) ds.c[0].pc[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
						if (!ds.c[1].pc[m]) ds.c[1].pc[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
						if (!ds.c[2].pc[m]) ds.c[2].pc[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
						if (run_cfg.nonlinear.seperate_linear_nonlinear)
						{
							if (!ds.c[0].pc_nl[m]) ds.c[0].pc_nl[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
							if (!ds.c[1].pc_nl[m]) ds.c[1].pc_nl[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); //subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
							if (!ds.c[2].pc_nl[m]) ds.c[2].pc_nl[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
							if (!ds.c[0].pc_prev_nl[m]) ds.c[0].pc_prev_nl[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
							if (!ds.c[1].pc_prev_nl[m]) ds.c[1].pc_prev_nl[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
							if (!ds.c[2].pc_prev_nl[m]) ds.c[2].pc_prev_nl[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
						}
						if (!ds.c[0].pc_prev[m]) ds.c[0].pc_prev[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
						if (!ds.c[1].pc_prev[m]) ds.c[1].pc_prev[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
						if (!ds.c[2].pc_prev[m]) ds.c[2].pc_prev[m] = fdtd::field4<double>(mp.nr_cp[m], subdomain.cells_x, subdomain.cells_y, subdomain.cells_z); // subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
					}
				}
			}
		}
		if (!subdomain.cpml_only)
		{
			// the volume DFT, and P_NL / DFG copies
			fdtd::monitors::volume_dft().allocate();
			if (run_cfg.thermal.save_output_temperature || run_cfg.thermal.load_temperature)
			{
				fdtd::materials::thermal::electron_temp_tf = fdtd::field3<double>(subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z);
			}
			fdtd::materials::material_map().allocate_inner(subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z, subdomain.cells_x, subdomain.cells_y, subdomain.cells_z);


		}

		return true;
	}
	catch (const std::bad_alloc &)
	{

		return false;
	}
}
/**
 * @brief Allocates the memory for the DFT fields and material index array which
 * can be used for 2D visualization
 *
 * @return true when memory allocation has been succesfull
 * @return false when error occurs
 */
bool fdtd::app::allocate_planes_and_material_maps()
{
	// the display DFT planes, on the animation planes of this rank, subdomain.plots_*()
	fdtd::monitors::dft_planes().allocate();
	return true;
}
/**
 * @brief Allocates the memory for the Signal which determines the injected
 * incident wave
 *
 * @return true when memory allocation has been succesfull
 * @return false when "bad_alloc" has been thrown
 */
bool fdtd::app::allocate_aux_source()
{
	fdtd::sources::AuxGrid &g = fdtd::sources::aux_grid();
	try
	{
		g.mat = fdtd::vector1<int>(run_cfg.grid.nr_cells_uni + 1);
		for (int j = 0; j <= run_cfg.grid.nr_cells_uni; ++j)
			g.mat[j] = j < run_cfg.grid.nr_cells_uni / 2 ? run_cfg.materials.first_medium : run_cfg.materials.second_medium;
		// the lines only on ranks with a TF/SF face
		if (subdomain.pw_nea_x || subdomain.pw_far_x || subdomain.pw_nea_y || subdomain.pw_far_y || subdomain.pw_nea_z || subdomain.pw_far_z)
		{
			g.line[0].allocate(run_cfg.grid.nr_cells_uni, run_cfg.grid.num_cpml_cells_y);
			if (run_cfg.grid.pbc_x && run_cfg.grid.pbc_z)
				g.line[1].allocate(run_cfg.grid.nr_cells_uni, run_cfg.grid.num_cpml_cells_y);
		}
	}
	catch (std::bad_alloc &)
	{
		return false;
	}
	return true;
}

/**
 * @brief Allocates the memory for the total field and scattered
 *  field monitors used for ART calculations
 *
 * @return true when memory allocation has been succesfull
 * @return false when error occurs
 */
bool fdtd::app::allocate_scattering_and_absorption()
{

	// the reflection and transmission planes of periodic runs
	fdtd::monitors::art_planes().allocate();
	// the DFT planes of the scattering and total flux boxes
	fdtd::monitors::flux_boxes().allocate();

	return true;
}
/**
 * @brief Allocates the dispersive state, polarisations, currents, critical-point poles, once the mesh
 *        is built, limiting the materials of boxed_state() to the rows (i, j) of their cells
 *
 * A material that can be a background, first, middle or second medium, or on the 1D layer profile
 * reaches the CPML cells and keeps full rows; the others occur only through the material maps of the
 * inner region, so their box is the bounding box of their cells there.
 * The rows outside a box read 0, as the full arrays did there, only a material's own cells are written.
 */
bool fdtd::app::allocate_material_state()
{
	auto &ds = fdtd::materials::dispersive_state();
	const auto &mp = fdtd::materials::model_parameters();
	const auto &S = subdomain;
	std::vector<bool> background(max_slots, false);
	auto mark = [&](int s) {
		if (s >= 0 && s < max_slots)
			background[s] = true;
	};
	mark(run_cfg.materials.first_medium);
	mark(run_cfg.materials.middle_medium);
	mark(run_cfg.materials.second_medium);
	for (int j = 0; j < fdtd::sources::nr_cell_uni + 1; ++j)
		mark(background_y[j]);

	// bounding boxes of the materials' inner cells, in cell indices
	std::vector<CellBox> box(max_slots, CellBox{1 << 30, -1, 1 << 30, -1});
	if (!S.cpml_only && S.inner_dim_z >= 0)
		box = fdtd::materials::material_map().dispersive_boxes(S.inner_org_x, S.inner_org_y);
	auto box_of = [&](int mm) -> std::optional<CellBox> {
		if (!boxed_state(mm) || background[mm])
			return std::nullopt;
		if (box[mm].i1 < 0)
			return CellBox{0, 0, 0, 0}; // not on this rank: one row
		return box[mm];
	};
	try
	{
		fdtd::materials::dispersive_state().for_each_field([&](auto &field) { field.commit(box_of); });
		for (int mm : deferred_poles)
		{
			const std::optional<CellBox> b = box_of(mm);
			const std::array<size_t, 4> shape = {size_t(mp.nr_cp[mm] > 0 ? mp.nr_cp[mm] : 1), size_t(S.cells_x + 1), size_t(S.cells_y + 1), size_t(S.cells_z + 1)};
			fdtd::materials::dispersive_state().for_each_pole(mm, [&](fdtd::array4<double> &pole) {
				pole = b ? fdtd::array4<double>::boxed(shape, b->i0, b->i1, b->j0, b->j1) : fdtd::field4<double>(mp.nr_cp[mm], S.cells_x, S.cells_y, S.cells_z);
			});
		}
		deferred_poles.clear();
	}
	catch (const std::bad_alloc &)
	{
		return false;
	}

	// memory of the state: stored, and what full rows would take
	double bytes[2] = {0, 0};
	for (auto *field : {&ds.c[0].p, &ds.c[0].p_prev, &ds.c[1].p, &ds.c[1].p_prev, &ds.c[2].p, &ds.c[2].p_prev})
		bytes[0] += double(field->stored()) * sizeof(double);
	for (int mm = 0; mm < max_slots; ++mm)
		fdtd::materials::dispersive_state().for_each_pole(mm, [&](fdtd::array4<double> &pole) { bytes[0] += double(pole.stored()) * sizeof(double); });
	const double cell_bytes = double(S.cells_x + 1) * (S.cells_y + 1) * (S.cells_z + 1) * sizeof(double);
	for (int mm = 0; mm < max_slots; ++mm)
		if (ds.c[0].p[mm])
			bytes[1] += cell_bytes * 6 * (1 + (ds.c[0].pc[mm] ? double(ds.c[0].pc[mm].extent(0)) : 0.0));
	double total[2] = {0, 0};
	MPI_Reduce(bytes, total, 2, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
	if (fdtd::grid::world_rank == ROOTPROCESS && total[1] > 0)
		fdtd::io::log("setup", "Dispersive state: %.1f MB (%.1f MB with full rows)\n", total[0] / 1048576.0, total[1] / 1048576.0);
	return true;
}

/// the rows outside the material boxes still read 0, a write there would be a bug; prints the arrays that fail
bool fdtd::app::material_state_intact()
{
	auto &ds = fdtd::materials::dispersive_state();
	bool ok = true;
	for (auto *field : {&ds.c[0].p, &ds.c[1].p_prev, &ds.c[0].p_prev, &ds.c[1].p, &ds.c[2].p, &ds.c[2].p_prev})
		ok = ok && field->padding_zero();
	for (int mm = 0; mm < max_slots; ++mm)
		fdtd::materials::dispersive_state().for_each_pole(mm, [&](fdtd::array4<double> &pole) { ok = ok && pole.padding_zero(); });
	if (!ok)
		std::fprintf(stderr, "rank %d: dispersive state written outside its material box\n", fdtd::grid::world_rank);
	return ok;
}

/// The cell-size reciprocals and the time step from the Courant factor and the smallest cell size.
void fdtd::app::set_time_step()
{
	fdtd::app::RunConfig &run_cfg = fdtd::app::edit_run_config(); // configure phase
	double ddd;
	if ((run_cfg.grid.dx != run_cfg.grid.dy) || (run_cfg.grid.dx != run_cfg.grid.dz) || (run_cfg.grid.dy != run_cfg.grid.dz))
	{
		if (run_cfg.grid.dx < run_cfg.grid.dy && run_cfg.grid.dx < run_cfg.grid.dz)
		{
			ddd = run_cfg.grid.dx;
		}
		else if (run_cfg.grid.dy < run_cfg.grid.dx && run_cfg.grid.dy < run_cfg.grid.dz)
		{
			ddd = run_cfg.grid.dy;
		}
		else
		{
			ddd = run_cfg.grid.dz;
		}
	}
	else
	{
		ddd = run_cfg.grid.dx = run_cfg.grid.dy = run_cfg.grid.dz;
	}

	run_cfg.time.dt = run_cfg.time.cfl_coeff * ddd / C0;
}
