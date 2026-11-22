/**
 * @file app/state_file.cc
 * @brief The state dump of an interrupted run and its restart, inc/app/state_file.h
 */
#include "app/state_file.h"
#include "app/run_config.h"
#include "grid/state.h"
#include "materials/state.h"
#include "fields/fields.h"
#include "boundary/cpml.h"
#include "app/archive.h"
#include "monitors/art_planes.h"
#include "monitors/custom_boxes.h"
#include "monitors/dft_planes.h"
#include "monitors/field_trackers.h"
#include "monitors/flux_boxes.h"
#include "monitors/state.h"
#include "monitors/time_probes.h"
#include "monitors/volume_dft.h"
#include "sources/aux_grid.h"
using fdtd::grid::subdomain;

namespace
{
const fdtd::app::RunConfig &run_cfg = fdtd::app::run_config(); ///< the run configuration

void plane(fdtd::io::Archive &ar, fdtd::monitors::DftPlane &p)
{
	ar.array(p.re);
	ar.array(p.im);
}

/// everything that changes during the time loop, in one fixed order
void transfer(fdtd::io::Archive &ar, int &step)
{
	// the run this file belongs to
	ar.same(std::uint64_t(0x3145544154534446ULL), "not a state file"); // "FDSTATE1"
	ar.same(fdtd::grid::world_size, "number of MPI ranks");
	ar.same(fdtd::grid::world_rank, "rank");
	ar.same(subdomain.cells_x, "subdomain cells x");
	ar.same(subdomain.cells_y, "subdomain cells y");
	ar.same(subdomain.cells_z, "subdomain cells z");
	ar.same(run_cfg.time.dt, "time step");
	ar.same(run_cfg.analysis.nr_an_frequencies, "number of analysed frequencies");
	ar.value(step);
	ar.value(fdtd::monitors::count_art);

	// fields, CPML, the 1D auxiliary grid
	auto &f = fdtd::fields::fields();
	ar.array(f.E);
	ar.array(f.H);
	ar.array(f.E_prev);
	ar.array(f.E_prev_2);
	ar.array(f.H_prev);
	ar.array(f.E_scratch);
	ar.array(fdtd::boundary::cpml().psi);
	for (auto &line : fdtd::sources::aux_grid().line)
	{
		ar.array(line.e);
		ar.array(line.e_prev);
		ar.array(line.h);
		ar.array(line.p);
		ar.array(line.p_prev);
		ar.array(line.pc);
		ar.array(line.pc_prev);
		ar.array(line.psi_e);
		ar.array(line.psi_h);
	}

	// dispersive state
	fdtd::materials::dispersive_state().for_each_saved([&](auto &x) { ar.array(x); });

	// monitors
	for (auto &box : fdtd::monitors::flux_boxes().box)
		for (auto &normal : box.face)
			for (auto &face : normal)
			{
				ar.array(face.re);
				ar.array(face.im);
			}
	auto &planes = fdtd::monitors::dft_planes();
	for (auto &p : planes.plane)
		plane(ar, p);
	plane(ar, planes.xz_1);
	plane(ar, planes.xz_2);
	plane(ar, planes.graphene);
	for (auto &p : fdtd::monitors::art_planes().plane)
		plane(ar, p);
	auto &volume = fdtd::monitors::volume_dft();
	ar.array(volume.field), ar.array(volume.pnl);
	ar.array(fdtd::monitors::field_trackers().dft);
	ar.array(fdtd::monitors::time_probes().E), ar.array(fdtd::monitors::time_probes().H);
	for (auto &box : fdtd::monitors::custom_boxes().frequency)
		ar.array(box.dft);
}

std::string rank_file(const std::string &dir) { return dir + "/rank_" + std::to_string(fdtd::grid::world_rank) + ".bin"; }
} // namespace

std::string fdtd::app::state_unsupported()
{
	if (run_cfg.materials.hydrodynamics)
		return "the hydrodynamic model";
	if (run_cfg.thermal.two_temperature_model)
		return "the two-temperature model";
	if (run_cfg.source.oblique_incidence_with_medium)
		return "the layered-medium wave source";
	if (run_cfg.source.mode_source)
		return "the mode source";
	return "";
}

bool fdtd::app::save_state(const std::string &dir, int step, std::string &why)
{
	try
	{
		fdtd::io::Archive ar(rank_file(dir), true);
		transfer(ar, step);
	}
	catch (const std::exception &e)
	{
		why = e.what();
		return false;
	}
	return true;
}

bool fdtd::app::load_state(const std::string &dir, int &step, std::string &why)
{
	try
	{
		fdtd::io::Archive ar(rank_file(dir), false);
		transfer(ar, step);
	}
	catch (const std::exception &e)
	{
		why = e.what();
		return false;
	}
	return true;
}
