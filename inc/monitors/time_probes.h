/**
 * @file monitors/time_probes.h
 * @brief The time-domain probes, [probe.<name>] sections: E and H at single cells, every step
 *
 * @verbatim
      point[p]          the probes inside this rank's output block: local index, with the inner
                        origin, and global index, in global (x, y, z) scan order
      E[c], H[c]        [probe][step], steps 0..N, 0 stays zero; E and H as they are after the E update
   @endverbatim
 * Written at every checkpoint to PROBES/TimePlot_<x>_<y>_<z>.txt: the material indices of the cell,
 * then one line per recorded step (Ex Ey Ez Hx Hy Hz). Recorded from the first step, also during
 * the settling periods of continuous-wave runs.
 */
#ifndef FDTD_MONITORS_TIME_PROBES_H
#define FDTD_MONITORS_TIME_PROBES_H

#include "core/array.h"
#include "monitors/inputs.h"
#include <vector>

namespace fdtd::monitors
{
struct ProbeCell
{
	int x; ///< local index
	int y; ///< local index
	int z; ///< local index
	int gx; ///< global index, as given in the INI
	int gy; ///< global index, as given in the INI
	int gz; ///< global index, as given in the INI
};

struct TimeProbes
{
	std::vector<ProbeCell> point;
	array2<double> E[3];
	array2<double> H[3];

	/// keeps the probes of this rank, none on ranks that are all CPML, and allocates their traces
	void setup(const std::vector<ProbePoint> &probes);
	void release();
	/// records E and H at step
	void sample(int step);
	/// writes the traces up to step - 1
	void write(int step) const;
};

TimeProbes &time_probes();
} // namespace fdtd::monitors

#endif
