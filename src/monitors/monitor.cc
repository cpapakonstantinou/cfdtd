/**
 * @file monitors/monitor.cc
 * @brief The monitor list, monitors/monitor.h
 */
#include "monitors/monitor.h"
#include "monitors/state.h"
#include "grid/subdomain.h"
#include "grid/state.h"
#include "sources/state.h"
#include "io/binary.h"
#include "grid/output_layout.h"
#include "monitors/art.h"
#include "monitors/art_planes.h"
#include "monitors/custom_boxes.h"
#include "monitors/dft_planes.h"
#include "monitors/exact_multipoles.h"
#include "monitors/field_points.h"
#include "monitors/field_trackers.h"
#include "monitors/flux_boxes.h"
#include "monitors/near_to_far.h"
#include "monitors/q_coefficients.h"
#include "monitors/time_probes.h"
#include "monitors/volume_dft.h"
#include <cmath>
#include <string>
#include "monitors/spec.h"
using fdtd::constants::C0;
using fdtd::constants::PI2;
using fdtd::grid::subdomain;

namespace
{
const fdtd::monitors::Spec &sp = fdtd::monitors::spec();        ///< what the monitors need from the run
} // namespace
using namespace fdtd::monitors; // transitional: module state, monitors/state.h

namespace fdtd::monitors
{
namespace
{
class VolumeFields : public Monitor
{
public:
	VolumeFields(bool inner, bool with_h) : _inner(inner), _with_h(with_h) {}
	void sample(const SampleTime &time) override
	{
		if (!_inner)
			return;
		if (_with_h)
			volume_dft().sample_h(time.h);
		volume_dft().sample_e(time.e);
	}
	void normalize() override
	{
		if (!_inner)
			return;
		volume_dft().normalize(0);
		if (_with_h)
			volume_dft().normalize(1);
	}
	/// collective; returns at once unless the exact multipoles are on
	void finish() override { write_exact_multipoles(); }
	/// the normalised volume E and H, one file per frequency, field saves
	void write_final() override
	{
		if (!sp.save_output || !_inner)
			return;
		grid::OutputLayout L = grid::output_layout();
		for (int m = 0; m < sp.nr_frequencies; m++)
		{
			const std::string lambda = std::to_string(std::lround((PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9)));
			save_as_binary_3d(sp.fields_dir + "/E_" + lambda + "nm_3D.h5", volume_dft().field[0][0][m], volume_dft().field[0][1][m], volume_dft().field[0][2][m], fdtd::grid::output_comm,
						   L.start, L.global, L.local, L.shift, false, true, false);
			save_as_binary_3d(sp.fields_dir + "/H_" + lambda + "nm_3D.h5", volume_dft().field[1][0][m], volume_dft().field[1][1][m], volume_dft().field[1][2][m], fdtd::grid::output_comm,
						   L.start, L.global, L.local, L.shift, false, true, false);
		}
	}

private:
	bool _inner;
	bool _with_h;
};

class NonlinearPolarisation : public Monitor
{
public:
	explicit NonlinearPolarisation(bool inner) : _inner(inner) {}
	void sample(const SampleTime &time) override
	{
		if (_inner)
			sample_nonlinear_polarisation(time.e);
	}
	void normalize() override
	{
		if (_inner)
			volume_dft().normalize_pnl();
	}

private:
	bool _inner;
};

class FluxBoxMonitor : public Monitor
{
public:
	explicit FluxBoxMonitor(bool inner) : _inner(inner) {}
	void sample(const SampleTime &time) override
	{
		if (_inner)
			flux_boxes().sample(time.e, time.h);
	}
	void normalize() override
	{
		if (_inner)
			flux_boxes().normalize();
	}
	void finish() override
	{
		write_q_coefficients();
		MPI_Barrier(MPI_COMM_WORLD);
		write_art_tfsf();
	}
	void checkpoint(int) override
	{
		write_q_coefficients();
		MPI_Barrier(MPI_COMM_WORLD);
		write_art_tfsf();
		MPI_Barrier(MPI_COMM_WORLD);
	}

private:
	bool _inner;
};

/// the scattering-box faces written for the near-to-far transform, SHG/THG runs with two frequencies
class NearToFar : public Monitor
{
public:
	void sample(const SampleTime &) override {}
	void finish() override { write_near_to_far(); }
};

/// the near-, far-field and forward/backward point files, read from the display and ART planes
class FieldPoints : public Monitor
{
public:
	void sample(const SampleTime &) override {}
	void finish() override { write_field_points(); }
};

class DisplayPlanes : public Monitor
{
public:
	void sample(const SampleTime &time) override { dft_planes().sample(time.e); }
	void normalize() override { dft_planes().normalize(); }
	/// the frequency-domain maps on the animation planes: with field saves the normalised volume E
	/// cut along the planes, otherwise the plane DFTs themselves
	void finish() override
	{
		if (!sp.dft_plot)
			return;
		if (!sp.save_output)
		{
			dft_planes().write_maps();
			return;
		}
		if (subdomain.cpml_only)
			return;
		grid::OutputLayout L = grid::output_layout();
		for (int m = 0; m < sp.nr_frequencies; m++)
		{
			const std::string lambda = std::to_string(std::lround((PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9)));
			const std::string base = sp.fields_dir + "/E_" + lambda;
			if (subdomain.plots_xy() && dft_plane_xy)
				save_as_binary_2d(base + "nm_XY.h5", volume_dft().field[0][0][m], volume_dft().field[0][1][m], volume_dft().field[0][2][m], fdtd::grid::output_comm_plot_xy, subdomain.k_animation,
							   L.start, L.global, L.local, L.shift, 0, 1, false, true, false);
			if (subdomain.plots_xz() && dft_plane_xz)
				save_as_binary_2d(base + "nm_XZ.h5", volume_dft().field[0][0][m], volume_dft().field[0][1][m], volume_dft().field[0][2][m], fdtd::grid::output_comm_plot_xz, subdomain.j_animation,
							   L.start, L.global, L.local, L.shift, 0, 2, false, true, false);
			if (subdomain.plots_yz() && dft_plane_yz)
				save_as_binary_2d(base + "nm_YZ.h5", volume_dft().field[0][0][m], volume_dft().field[0][1][m], volume_dft().field[0][2][m], fdtd::grid::output_comm_plot_yz, subdomain.i_animation,
							   L.start, L.global, L.local, L.shift, 1, 2, false, true, false);
		}
	}
};

class ArtPlaneMonitor : public Monitor
{
public:
	void sample(const SampleTime &time) override { art_planes().sample(time.e, time.h); }
	void normalize() override { art_planes().normalize(); }
	void finish() override
	{
		if (sp.art)
			write_art_spectra(_last_step);
	}
	/// the intermediate spectra at steps N - kN/10, k = 9..1, count_art numbers them
	void after_step(int step) override
	{
		_last_step = step;
		if (sp.grating || !sp.art)
			return;
		for (int k = 1; k <= 9; ++k)
			if (step == sp.steps - (k * sp.steps) / 10)
			{
				write_art_spectra(step);
				count_art += 1;
				return;
			}
	}

private:
	int _last_step = 0; ///< the last step after_step saw: the step the final spectra are taken at
};

class FieldTrackerMonitor : public Monitor
{
public:
	void sample(const SampleTime &time) override { field_trackers().sample(time.e, time.h); }
	void checkpoint(int) override { field_trackers().write(); }
};

class CustomMonitors : public Monitor
{
public:
	void sample(const SampleTime &time) override { custom_boxes().sample(time.e, time.step); }
	void write_final() override { custom_boxes().write_dft(); }
};

class TimeProbeMonitor : public Monitor
{
public:
	void sample(const SampleTime &time) override { time_probes().sample(time.step); }
	bool every_step() const override { return true; }
	void checkpoint(int step) override { time_probes().write(step); }
};
} // namespace

std::vector<std::unique_ptr<Monitor>> build_monitors()
{
	// membership from run-wide flags only, finish is collective; the finish order is the old one:
	// exact multipoles, ART, near-to-far, Q coefficients, field points
	std::vector<std::unique_ptr<Monitor>> list;
	const bool inner = !subdomain.cpml_only;
	if (sp.save_output || sp.abs_sigma || multipoles_wanted())
		list.push_back(std::make_unique<VolumeFields>(inner, sp.save_output));
	if (multipoles_wanted() && (sp.shg || sp.thg))
		list.push_back(std::make_unique<NonlinearPolarisation>(inner));
	list.push_back(std::make_unique<ArtPlaneMonitor>());
	if ((sp.shg || sp.thg) && sp.nr_frequencies == 2)
		list.push_back(std::make_unique<NearToFar>());
	if (q_coeff)
		list.push_back(std::make_unique<FluxBoxMonitor>(inner));
	list.push_back(std::make_unique<FieldPoints>());
	list.push_back(std::make_unique<DisplayPlanes>());
	if (sp.sphere_multipoles) // the tracking points, some ranks hold none
		list.push_back(std::make_unique<FieldTrackerMonitor>());
	if (custom_boxes().any()) // built from the configuration: the same on every rank
		list.push_back(std::make_unique<CustomMonitors>());
	list.push_back(std::make_unique<TimeProbeMonitor>()); // no-op on ranks without probes
	return list;
}
} // namespace fdtd::monitors
