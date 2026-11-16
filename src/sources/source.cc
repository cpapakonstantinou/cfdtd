/**
 * @file sources/source.cc
 * @brief The list of this rank's sources, sources/source.h
 */
#include "sources/source.h"
#include "sources/state.h"
#include "sources/signal.h"
#include "sources/one_dim_simulator.h"
#include "core/views.h"
#include "grid/subdomain.h"
#include "fields/fields.h"
#include "grid/state.h"
#include "sources/oblique.h"
#include "sources/superposition.h"
#include "sources/tfsf.h"
#include "sources/tfsf_source.h"
#include "sources/mode_source.h"
#include "grid/lattice.h"
#include "sources/spec.h"
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::fields::Hx;
using fdtd::fields::Hy;
using fdtd::fields::Hz;
using fdtd::grid::subdomain;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::sources; // transitional: module state, sources/state.h

namespace fdtd::sources
{
namespace
{
/// mode source: the mode profile times the signal on the near y face, sources/tfsf.h
class ModeSource : public Source
{
public:
	explicit ModeSource(const Signal &signal) : _signal(signal) {}
	void inject_h(const StepTime &time) const override { mode_source_h(_signal, time.t); }
	void inject_e(const StepTime &time) const override { mode_source_e(_signal, time.t); }

private:
	const Signal &_signal;
};

/// analytic plane-wave superposition or Gaussian beam, sources/superposition.h
class Superposition : public Source
{
public:
	Superposition(const Signal &signal, const std::vector<PlaneWave> &waves, const WaveFrame &frame) : _signal(signal), _waves(waves), _frame(frame) {}
	void inject_h(const StepTime &time) const override { superposition_h(_signal, _frame, time.h, _waves); }
	void inject_e(const StepTime &time) const override { superposition_e(_signal, _frame, time.e, _waves); }

private:
	const Signal &_signal;
	const std::vector<PlaneWave> &_waves;
	WaveFrame _frame; ///< focus, TF/SF faces, common delay
};

/// plane waves through a layered background, from 1D simulations, sources/oblique.h
class LayeredSuperposition : public Source
{
public:
	explicit LayeredSuperposition(OneDimensionalSimulator &ods) : _ods(ods) {}
	void inject_h(const StepTime &time) const override { oblique_h(_ods, time.h); }
	void inject_e(const StepTime &time) const override { oblique_e(_ods, time.e); }

private:
	OneDimensionalSimulator &_ods;
};

/**
 * @brief The point dipoles that lie on this rank
 * @verbatim
      E_c[cell] += signal(step dt) * e_c        (c = x, y, z)
      H_c[cell] += signal(step dt) * 0, no magnetic dipoles in the configuration; the zero
                                                  term is kept so signed zeros match the previous build
   @endverbatim
 * A dipole on a plane shared by two ranks is injected by both. The signal is tabulated per step when
 * the source is built, before: a [step][dipole][6] table of products.
 */
class Dipoles : public Source
{
public:
	Dipoles(const Signal &signal, int steps)
	{
		for (const auto &d : dipole_input)
			if (d.x >= subdomain.output_start[0] && d.x < subdomain.output_start[0] + subdomain.output_size[0] && d.y >= subdomain.output_start[1] && d.y < subdomain.output_start[1] + subdomain.output_size[1] &&
				d.z >= subdomain.output_start[2] && d.z < subdomain.output_start[2] + subdomain.output_size[2])
				_dipoles.push_back({int(d.x - subdomain.output_start[0] + subdomain.inner_org_x), int(d.y - subdomain.output_start[1] + subdomain.inner_org_y),
									int(d.z - subdomain.output_start[2] + subdomain.inner_org_z), {d.e_x, d.e_y, d.e_z}});
		if (!_dipoles.empty())
		{
			_signal.resize(steps + 1); // steps 1 .. N, the previous table had N rows: step N read past it
			for (int j = 0; j <= steps; j++)
				_signal[j] = signal.value(j * sp.dt);
		}
	}
	void inject_h(const StepTime &time) const override
	{
		const std::array<double, 3> none{0, 0, 0};
		for (const Dipole &d : _dipoles)
			add(d, {Hx, Hy, Hz}, none, time.step);
	}
	void inject_e(const StepTime &time) const override
	{
		for (const Dipole &d : _dipoles)
			add(d, {Ex, Ey, Ez}, d.e, time.step);
	}

private:
	struct Dipole
	{
		int i; ///< local cell
		int j; ///< local cell
		int k; ///< local cell
		std::array<double, 3> e;   ///< E amplitudes
	};
	void add(const Dipole &d, const std::array<view3<double>, 3> &f, const std::array<double, 3> &amp, int step) const
	{
		for (int c = 0; c < 3; ++c)
			f[c][d.i][d.j][d.k] += _signal[step] * amp[c];
	}
	std::vector<Dipole> _dipoles;
	std::vector<double> _signal; ///< signal at each step
};
} // namespace

SourceList build_sources(const Signal &signal, const std::vector<PlaneWave> &waves, const WaveFrame &frame, OneDimensionalSimulator &ods, int steps)
{
	SourceList list;
	const bool superposition = subdomain.has_plane_wave && (sp.gauss_beam || sp.superposition) && !lat.periodic_x && !lat.periodic_z;
	if (subdomain.has_plane_wave && sp.tfsf && !sp.mode_source)
		list.push_back(std::make_unique<TfsfSource>(TfsfSource::for_this_rank(signal)));
	if (sp.mode_source)
		list.push_back(std::make_unique<ModeSource>(signal));
	if (superposition && !sp.oblique_with_medium)
		list.push_back(std::make_unique<Superposition>(signal, waves, frame));
	if (superposition && sp.oblique_with_medium)
		list.push_back(std::make_unique<LayeredSuperposition>(ods));
	if (sp.dipole)
		list.push_back(std::make_unique<Dipoles>(signal, steps));
	return list;
}
} // namespace fdtd::sources
