#ifndef __DRIVER_CONTEXT_H__
#define __DRIVER_CONTEXT_H__
/** \file   app/driver_context.h
 *  \brief  Private, legacy-typed state of fdtd::Driver, shared by driver.cc and time_loop.cc.
 **/
#include "grid/output_layout.h"
#include "monitors/monitor.h"
#include "sources/source.h"
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "core/array.h"
#include "app/driver.h"

#include "io/binary.h"
#include "sources/one_dim_simulator.h"

namespace fdtd
{

struct Driver::Context
{
	// subdomain bookkeeping for the binary writers
	grid::OutputLayout layout{};

	std::vector<PlaneWave> waves; ///< plane waves of the superposition / Gaussian beam, sources/waves.h
	WaveFrame frame;               ///< where they are anchored, focus, TF/SF faces, delay
	OneDimensionalSimulator ods;   ///< their 1D simulations through a layered background

	std::string movie_base;   ///< <output>/frames/Efield_, then <plane>_<frame number>.h5

	// per sweep
	/// time-domain probe recordings [probe][timestep]; released automatically, were leaked before

	/// this rank's sources, built once per sweep, sources/source.h
	fdtd::sources::SourceList sources;
	std::vector<std::unique_ptr<fdtd::monitors::Monitor>> monitors;

	bool   output_pump_data = false;
	double loop_start_wtime = 0.0;
	int    stop_at_step = -1;      ///< wall-clock limit, -1 = none
	int    checkpoint_period = 1;  ///< max(1, timesteps / [run] checkpoint_divisions)
};

} // namespace fdtd

#endif // __DRIVER_CONTEXT_H__
