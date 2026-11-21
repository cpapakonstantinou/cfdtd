#include "sources/types.h"
#include "io/text.h"
#include "grid/subdomain.h"
#include "core/constants.h"
#include "materials/state.h"
#include "grid/state.h"
#include "sources/one_dim_simulator.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <string>
#include <iostream>
#include "grid/lattice.h"
#include "sources/spec.h"
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace

bool OneDimensionalSimulator::init(std::vector<PlaneWave> planeWaves, const WaveFrame &frame)
{
	_frame = frame;
	bool taskHasTfsfSurface = subdomain.pw_nea_x || subdomain.pw_nea_y || subdomain.pw_nea_z || subdomain.pw_far_x || subdomain.pw_far_y || subdomain.pw_far_z;
	bool taskRunsOneDimSimulation = taskHasTfsfSurface || (fdtd::grid::world_rank == ROOTPROCESS);

	if (!taskRunsOneDimSimulation)
		return false;

	check_media_ordering();

	for (PlaneWave &tmpWave : planeWaves)
	{
		PlaneWaveOneDimensionalSimulator singleSimulator;
		singleSimulator.init(tmpWave, frame);
		singleSimulatorList.push_back(singleSimulator);
	}

	if (fdtd::grid::world_rank==ROOTPROCESS)
	{
		write_onedim_results();
	}

	return true;
}

bool OneDimensionalSimulator::write_onedim_results()
{
	static const char *name[6] = {"ex", "ey", "ez", "hx", "hy", "hz"};
	FILE *f[6];
	for (int c = 0; c < 6; ++c)
		f[c] = fdtd::io::open_report(sp.reports_dir + "/aux_grid_" + name[c] + ".txt",
									 "one row per frame period, one column per cell of the auxiliary grid along y");

	std::array<double, 6> fields;

	for (int i = 0; i < sp.one_dim_steps; i++)
	{
		int plotSteps = std::max({sp.plot_yz, sp.plot_xy, sp.plot_xz});
		if (i % plotSteps == 0)
		{
			double t = i * sp.dt;

			int gridShift = singleSimulatorList[0].gridShiftIdx;
			int gridLength = singleSimulatorList[0].gridLength;

			// the last point is left out: its staggered interpolation would read one past the arrays
			for (int j = -gridShift; j < gridLength - gridShift - 1; j++)
			{
				double y = j * lat.d[1];
				fields = return_fields(0.0, y, 0.0, t);
				for (int c = 0; c < 6; ++c)
					if (f[c])
						std::fprintf(f[c], "%s ", std::to_string(fields[c]).c_str());
			}
			for (int c = 0; c < 6; ++c)
				if (f[c])
					std::fputc('\n', f[c]);
		}
	}
	for (int c = 0; c < 6; ++c)
		if (f[c])
			std::fclose(f[c]);

	return true;
}

std::array<double, 6> OneDimensionalSimulator::return_fields(double x, double y, double z,
															 double t)
{
	// initialize an emtpy fields array
	std::array<double, 6> fields = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	// iterate through all simulators for individual plane waves
	for (PlaneWaveOneDimensionalSimulator &singleSimulator : singleSimulatorList)
	{
		// get fields from single plane wave
		std::array<double, 6> tmpFields = singleSimulator.return_fields(x, y, z, t);

		// add fields
		for (int i = 0; i < int(fields.size()); i++)
			fields[i] += tmpFields[i];
	}
	return fields;
}

int OneDimensionalSimulator::get_medium_index(double y)
{
	// use the first planeWave ODS, since the structure is the same in every
	return singleSimulatorList[0].get_medium_index(y);
}

/**
 * @brief Check if media are ordered by increasing density medium.
 * Currently we only support going to media with higher refractive indices.
 */
void OneDimensionalSimulator::check_media_ordering()
{
	const auto &mp = fdtd::materials::model_parameters();
	if (mp.n_refr_index[sp.medium[1]] < mp.n_refr_index[sp.medium[0]] && sp.middle_layer > 0.0)
		fdtd::io::log("sources", "WARNING: Middle medium has lower refractive index than first medium. "
			   "This functionality is currently not implemented and might cause problems! \n");

	if (mp.n_refr_index[sp.medium[2]] < mp.n_refr_index[sp.medium[0]])
		fdtd::io::log("sources", "WARNING: Second medium has lower refractive index than first medium. "
			   "This functionality is currently not implemented and might cause problems! \n");

	if (mp.n_refr_index[sp.medium[2]] < mp.n_refr_index[sp.medium[1]] && sp.middle_layer > 0.0)
		fdtd::io::log("sources", "WARNING: Second medium has lower refractive index than middle medium. "
			   "This functionality is currently not implemented and might cause problems! \n");
}