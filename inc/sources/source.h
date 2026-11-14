/**
 * @file sources/source.h
 * @brief The sources of a run as one list of objects with injection hooks
 *
 * Each step the time loop calls inject_h after the H update and inject_e after the E update, on
 * every source of this rank in list order. The list is built once per sweep (build_sources) from
 * the run's configuration, so no flag is tested inside the loop.
 * @verbatim
      source                 H step                               E step
      TF/SF plane wave       signal into the 1D grid, H faces     E faces
      mode source            MPB profile, H on its plane          same, E, any axis
      superposition / beam   analytic waves on the six faces      same, E
      layered-medium waves   1D simulations on the six faces      same, E
      dipoles                precomputed H terms at the dipoles   E terms
   @endverbatim
 */
#ifndef FDTD_SOURCES_SOURCE_H
#define FDTD_SOURCES_SOURCE_H

#include <memory>
#include <vector>

class OneDimensionalSimulator;
struct PlaneWave;
struct WaveFrame;

namespace fdtd::sources
{
/// the times a source is evaluated at in one step
struct StepTime
{
	double t;  ///< time of the step, the time loop's t
	double h;  ///< injection time of the H step
	double e;  ///< injection time of the E step
	int step;  ///< time step, precomputed sources
};

class Source
{
public:
	virtual ~Source() = default;
	virtual void inject_h(const StepTime &time) const = 0;
	virtual void inject_e(const StepTime &time) const = 0;
};

using SourceList = std::vector<std::unique_ptr<const Source>>;

/// the sources of this rank, in injection order; waves: the superposition / Gaussian beam waves,
/// ods: their 1D simulations through a layered background, both kept by the caller
class Signal;
/// the sources of this rank: each holds the run's signal; the superposition sources also the waves and their frame
SourceList build_sources(const Signal &signal, const std::vector<PlaneWave> &waves, const WaveFrame &frame, OneDimensionalSimulator &ods,
						 int steps); ///< steps: the run's step count, the dipole signal table
} // namespace fdtd::sources

#endif
