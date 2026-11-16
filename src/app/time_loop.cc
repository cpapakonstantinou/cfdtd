/** \file   app/time_loop.cc
 *  \brief  The FDTD time loop of one sweep.
 *
 *  The former main() re-tested ~20 source flags inside every time step. The
 *  source configuration cannot change during a sweep, so build_run_lists()
 *  resolves it once into a list of source objects, sources/source.h, and the
 *  loop calls their inject_h / inject_e hooks with the step's times. The
 *  monitors are resolved the same way, monitors/monitor.h, and sampled once
 *  per step.
 **/
#include "driver_context.h"
#include "app/state_file.h"
#include "io/text.h"
#include "grid/state.h"
#include "materials/thermal/thermal.h"
#include "materials/hydro/hydro.h"
#include "sources/signal.h"
#include "update/update.h"
#include "materials/nonlinear/nonlinear.h"
#include "app/signals.h"
#include "grid/halo_exchange.h"
#include "fields/halos.h"
#include "monitors/field_trackers.h"
#include "monitors/monitor.h"
#include "sources/source.h"

#include <cmath>
#include <cstdio>
#include "grid/subdomain.h"
#include "app/run_config.h"
#include "app/run_state.h"
#include "materials/hydro/state.h"
#include "fields/fields.h"
#include "materials/thermal/state.h"
#include "monitors/state.h"
#include "geometry/state.h"
#include "sources/aux_grid.h"
using fdtd::constants::ROOTPROCESS;
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::fields::Hx;
using fdtd::fields::Hy;
using fdtd::fields::Hz;
using fdtd::fields::ex_nl;
using fdtd::fields::ey_nl;
using fdtd::fields::ez_nl;
using fdtd::fields::hx_nl;
using fdtd::fields::hy_nl;
using fdtd::fields::hz_nl;
using fdtd::grid::subdomain;

namespace
{
const fdtd::app::RunConfig &run_cfg = fdtd::app::run_config(); ///< the run configuration
} // namespace

// File-local state, formerly globals used only in this file.
namespace
{
double t_pw_for_h;                                            ///< Source time used for the H plane-wave correction this step in s.
double t_pw_for_e;                                            ///< Source time used for the E plane-wave correction this step in s.
int max_time_itr;                                           ///< Last time step of the loop, copy of nr_max_time_steps.
} // namespace

namespace fdtd
{

// Source injection

void Driver::build_run_lists()
{
	_ctx->sources = fdtd::sources::build_sources(fdtd::sources::run_signal(), _ctx->waves, _ctx->frame, _ctx->ods, fdtd::app::run_state().nr_max_time_steps);
	_ctx->monitors = fdtd::monitors::build_monitors();
}

// Field steps

void Driver::update_magnetic()
{
	// the 1D plane-wave grid of the TF/SF source steps first, the TF/SF terms read it
	if (subdomain.has_plane_wave && run_cfg.source.tfsf)
		fdtd::sources::aux_grid().update_h();
	fdtd::update::update_h();
	const fdtd::sources::StepTime time{fdtd::app::run_state().t, t_pw_for_h, t_pw_for_e, fdtd::app::run_state().timestep};
	for (const auto& source : _ctx->sources)
		source->inject_h(time);

	fields::exchange_halos(grid::Field::H, 0);
	if (run_cfg.nonlinear.seperate_linear_nonlinear)
	{
		fields::exchange_halos(grid::Field::H, 1);
	}

	if (run_cfg.thermal.two_temperature_model)
	{
		const bool pump_ready = fdtd::materials::thermal::update_ttm(fdtd::app::run_state().timestep);
		_ctx->output_pump_data = _ctx->output_pump_data || pump_ready;
	}
	if (run_cfg.materials.hydrodynamics)
		fdtd::materials::hydro::update_hydrodynamics();
}

void Driver::update_electric()
{
	auto& c = *_ctx;
	if (subdomain.has_plane_wave && run_cfg.source.tfsf)
		fdtd::sources::aux_grid().update_e();
	fdtd::update::update_e();

	if (!run_cfg.materials.hydrodynamics && run_cfg.analysis.multipoles)
		fdtd::materials::hydro::update_electron_density_local();
	if (run_cfg.materials.hydrodynamics && !run_cfg.thermal.two_temperature_model)
		fdtd::materials::hydro::update_electron_density();

	const fdtd::sources::StepTime time{fdtd::app::run_state().t, t_pw_for_h, t_pw_for_e, fdtd::app::run_state().timestep};
	for (const auto& source : c.sources)
		source->inject_e(time);

	// nonlinearity added to the E update: also exchange the longitudinal components
	if ((run_cfg.nonlinear.second_harmonic || run_cfg.nonlinear.third_harmonic) && !run_cfg.nonlinear.pnl_approach)
	{
		if (fdtd::grid::world_rank == ROOTPROCESS && fdtd::app::run_state().timestep == 1)
			fdtd::io::log("run", "Nonlinearity added to the electric-field update\n");
		fdtd::materials::nonlinear::update_field();
		fields::exchange_longitudinal_e();
	}

	fields::exchange_halos(grid::Field::E, 0);

	if (run_cfg.nonlinear.seperate_linear_nonlinear)
	{
		fields::exchange_halos(grid::Field::E, 1);
	}

	if (run_cfg.materials.nonlocal && !run_cfg.materials.hydrodynamics)
	{
		fdtd::grid::exchange_margins(fdtd::materials::hydro::px_nl, fdtd::materials::hydro::py_nl, fdtd::materials::hydro::pz_nl);
		fdtd::materials::hydro::update_nonlocal();
	}

}

// Monitors and output

void Driver::accumulate_monitors()
{
	// CW signals: the DFTs accumulate only after the settling periods
	const bool settled = !(run_cfg.source.continuous_wave && fdtd::app::run_state().t < run_cfg.source.start_cw * double(1 / run_cfg.source.fundamental_frequency));
	const fdtd::monitors::SampleTime time{fdtd::app::run_state().t, fdtd::app::run_state().t_1, fdtd::app::run_state().timestep};
	for (const auto& monitor : _ctx->monitors)
		if (settled || monitor->every_step())
			monitor->sample(time);
}

void Driver::print_trace()
{
	const bool point = subdomain.plots_xy() && subdomain.plots_yz() && subdomain.plots_xz();
	if (!point)
		return;
	const int I = subdomain.i_animation;
	const int J = subdomain.j_animation;
	const int K = subdomain.k_animation;

	if (_cfg.run.trace_interval <= 0 || fdtd::app::run_state().timestep % _cfg.run.trace_interval != 0)
		return;
	fdtd::io::log("trace", "step %d signal %.5e Ex %.5e Ey %.5e Ez %.5e", fdtd::app::run_state().timestep, fdtd::sources::signal_value(fdtd::app::run_state().t), Ex[I][J][K], Ey[I][J][K], Ez[I][J][K]);
	if (run_cfg.thermal.two_temperature_model)
		fdtd::io::log("trace", "Te %.5e Tl %.5e", fdtd::materials::thermal::electron_temp[I][J][K], fdtd::materials::thermal::lattice_temp[I][J][K]);
}

/// Collective wall-clock guard, was: hardcoded 23.5 h measured with CPU clock().
bool Driver::wallclock_exceeded()
{
	auto& c = *_ctx;
	if (_cfg.run.wallclock_limit_hours <= 0)
		return false;
	if (fdtd::app::run_state().timestep == _cfg.run.wallclock_check_step)
	{
		if (fdtd::grid::world_rank == ROOTPROCESS)
		{
			const double elapsed = MPI_Wtime() - c.loop_start_wtime;
			c.stop_at_step = int(std::ceil(_cfg.run.wallclock_limit_hours * 3600.0 * fdtd::app::run_state().timestep / std::max(elapsed, 1e-9)));
			fdtd::io::log("run", "wall-clock limit %.2f h -> stopping at step %d at the latest\n", _cfg.run.wallclock_limit_hours, c.stop_at_step);
		}
		MPI_Bcast(&c.stop_at_step, 1, MPI_INT, ROOTPROCESS, MPI_COMM_WORLD);
	}
	if (c.stop_at_step > 0 && fdtd::app::run_state().timestep >= c.stop_at_step)
	{
		if (fdtd::grid::world_rank == ROOTPROCESS)
			fdtd::io::log("run", "run interrupted at step %d (wall-clock limit)\n", fdtd::app::run_state().timestep);
		return true;
	}
	return false;
}

void Driver::write_time_domain_frames()
{
	auto& c = *_ctx;
	if (subdomain.cpml_only)
		return;
	auto frame = [](int period) { return period != 0 && fdtd::app::run_state().timestep % period == 0; };
	// <frames>/Efield_<plane>_<frame number>.h5: the number last, so viewers load the files as a series
	auto frame_name = [&](const char *plane, int period) {
		char number[16];
		std::snprintf(number, sizeof number, "%05d", fdtd::app::run_state().timestep / period);
		return c.movie_base + plane + "_" + number + ".h5";
	};

	if (subdomain.plots_xy() && frame(run_cfg.output.plot_xy))
		save_as_binary_2d(frame_name("XY", run_cfg.output.plot_xy), Ex, Ey, Ez, fdtd::grid::output_comm_plot_xy, subdomain.k_animation,
					   c.layout.start, c.layout.global, c.layout.local, c.layout.shift, 0, 1, false, true, true);
	if (subdomain.plots_yz() && frame(run_cfg.output.plot_yz))
		save_as_binary_2d(frame_name("YZ", run_cfg.output.plot_yz), Ex, Ey, Ez, fdtd::grid::output_comm_plot_yz, subdomain.i_animation,
					   c.layout.start, c.layout.global, c.layout.local, c.layout.shift, 1, 2, false, true, true);
	if (subdomain.plots_xz() && frame(run_cfg.output.plot_xz))
		save_as_binary_2d(frame_name("XZ", run_cfg.output.plot_xz), Ex, Ey, Ez, fdtd::grid::output_comm_plot_xz, subdomain.j_animation,
					   c.layout.start, c.layout.global, c.layout.local, c.layout.shift, 0, 2, false, true, true);
	if (frame(run_cfg.output.time_plot_3d))
		save_as_binary_3d(frame_name("3D", run_cfg.output.time_plot_3d), Ex, Ey, Ez, fdtd::grid::output_comm,
					   c.layout.start, c.layout.global, c.layout.local, c.layout.shift, false, true, true);
}

/// Work done a fixed number of times per run, plus the ART checkpoints.
void Driver::periodic_analysis()
{
	auto& c = *_ctx;
	int i;
	int j;
	int k;

	const int frame_period = std::max({run_cfg.output.plot_xy, run_cfg.output.plot_yz, run_cfg.output.plot_xz, run_cfg.output.time_plot_3d});
	if (run_cfg.thermal.save_output_temperature && frame_period != 0 && fdtd::app::run_state().timestep % frame_period == 0)
	{
		for (i = 0; i <= subdomain.inner_dim_x; ++i)
			for (j = 0; j <= subdomain.inner_dim_y; ++j)
				for (k = 0; k <= subdomain.inner_dim_z; ++k)
					fdtd::materials::thermal::electron_temp_tf[i][j][k] = fdtd::materials::thermal::electron_temp[i + subdomain.inner_org_x][j + subdomain.inner_org_y][k + subdomain.inner_org_z];
		fdtd::materials::thermal::save_temperature(fdtd::app::run_state().timestep);
	}

	// the monitors' own schedules, the intermediate ART spectra
	for (const auto& monitor : c.monitors)
		monitor->after_step(fdtd::app::run_state().timestep);

	if (fdtd::app::run_state().timestep % c.checkpoint_period != 0)
		return;

	// the monitors' checkpoint output, time-domain probes, field trackers, Q coefficients
	for (const auto& monitor : c.monitors)
		monitor->checkpoint(fdtd::app::run_state().timestep);

}

// The loop

void Driver::run_time_loop()
{
	auto& c = *_ctx;
	fdtd::app::run_state().prgstart = clock();
	c.loop_start_wtime = MPI_Wtime();
	c.stop_at_step = -1;
	c.output_pump_data = false;
	c.checkpoint_period = std::max(1, fdtd::app::run_state().nr_max_time_steps / _cfg.run.checkpoint_count); // legacy N/10 divided by zero for N < 10
	if (_completed_steps == 0)
		fdtd::monitors::count_art = 0;
	max_time_itr = fdtd::app::run_state().nr_max_time_steps;

	if (fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("run", "Time loop: %d steps, dt = %e s\n", max_time_itr, run_cfg.time.dt);

	bool stopped = false;
	for (fdtd::app::run_state().timestep = _completed_steps + 1; fdtd::app::run_state().timestep <= max_time_itr; ++fdtd::app::run_state().timestep)
	{
		fdtd::app::run_state().t = fdtd::app::run_state().timestep * run_cfg.time.dt;
		fdtd::app::run_state().t_1 = fdtd::app::run_state().t + double(0.5) * run_cfg.time.dt;
		t_pw_for_h = fdtd::app::run_state().t - run_cfg.time.dt;
		t_pw_for_e = fdtd::app::run_state().t - double(0.5) * run_cfg.time.dt;

		update_magnetic();
		update_electric();
		accumulate_monitors();
		print_trace();
		write_time_domain_frames();
		MPI_Barrier(MPI_COMM_WORLD);
		periodic_analysis();
		MPI_Barrier(MPI_COMM_WORLD);
		if ((_signal = app::agree_signal(MPI_COMM_WORLD)) != 0)
		{
			if (fdtd::grid::world_rank == ROOTPROCESS)
				fdtd::io::log("run", "Signal %d received: stopping after step %d of %d, writing the outputs\n", _signal, fdtd::app::run_state().timestep, max_time_itr);
			stopped = true;
		}
		if (wallclock_exceeded()) // the step is complete either way: a restart continues with the next one
			stopped = true;
		if (stopped)
			break;
	}
	if (stopped && fdtd::app::run_state().timestep < max_time_itr)
		save_state();

	agree(true, "");
}

/// Every rank writes its state file, so `simulate --restart` can continue the run.
void Driver::save_state()
{
	const std::string reason = fdtd::app::state_unsupported();
	if (!reason.empty())
	{
		if (fdtd::grid::world_rank == ROOTPROCESS)
			fdtd::io::log("state", "not written: restart is not available with %s", reason.c_str());
		return;
	}
	std::string why;
	const bool written = fdtd::app::save_state(run_cfg.output.state_dir, fdtd::app::run_state().timestep, why);
	int all = written;
	int every = 0;
	MPI_Allreduce(&all, &every, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	_state_written = every != 0;
	if (!written)
		fdtd::io::log("state", "not written: %s", why.c_str());
	else if (fdtd::grid::world_rank == ROOTPROCESS && every)
		fdtd::io::log("state", "written after step %d to %s: continue with simulate --restart", fdtd::app::run_state().timestep, run_cfg.output.state_dir.c_str());
}

/// Restart: the run is set up; load what the interrupted run left and continue after its last step.
void Driver::load_state()
{
	std::string why = fdtd::app::state_unsupported();
	int step = 0;
	bool loaded = false;
	if (!why.empty())
		why = "restart is not available with " + why;
	else
		loaded = fdtd::app::load_state(run_cfg.output.state_dir, step, why);
	if (!agree(loaded, why.empty() ? "restart failed on another rank" : why))
		return;
	_completed_steps = step;
	if (fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("state", "restart after step %d of %d", step, fdtd::app::run_state().nr_max_time_steps);
}

} // namespace fdtd
