/**
 * @file monitors/custom_boxes.h
 * @brief The custom monitor boxes, [monitor.<name>] sections: a box of cells recorded in frequency
 *        running DFT of E, or in time, E every `step` steps
 *
 * @verbatim
      global_min/max[a]   the box in the inner domain: ceil(center -/+ size/2), max exclusive
      local_min/max[a]    its part on this rank, inner indices, here = the rank holds a part
      relative_start[a]   offset of this rank's part in the box, the writers' file offset
      comm                the ranks holding a part, MPI_COMM_NULL elsewhere
      dft[c]              frequency boxes: [frequency][i][j][k] over the local part
   @endverbatim
 * Outputs:
 * @verbatim
      frequency box   DFT/E_DFT_<name>_WL_<nm>nm.h5 at the end of the run, DFT/E_DFT_<name>_Mask.h5
      time box        PROBES/E_Probe_<name>_Frame_<n>.h5 every `step` steps, PROBES/E_Probe_<name>_Mask.h5
   @endverbatim
 */
#ifndef FDTD_MONITORS_CUSTOM_BOXES_H
#define FDTD_MONITORS_CUSTOM_BOXES_H

#include "core/carray.h"
#include "mpi.h"
#include "monitors/inputs.h"
#include <complex>
#include <string>
#include <vector>

namespace fdtd::monitors
{
struct CustomBox
{
	std::string name;
	int step = 1; ///< time boxes: frame period
	int global_min[3] = {};
	int global_max[3] = {};
	int global_count[3] = {};
	bool here = false;
	int local_min[3] = {};
	int local_max[3] = {};
	int local_count[3] = {};
	int relative_start[3] = {};
	MPI_Comm comm = MPI_COMM_NULL;
	ctetrad<std::complex<double>> dft[3];
};

struct CustomBoxes
{
	std::vector<CustomBox> frequency; ///< in INI order
	std::vector<CustomBox> time; ///< in INI order

	/// builds the boxes of the configuration, collective: one communicator per box
	void build(const std::vector<MonitorBox> &boxes);
	void release();
	bool any() const { return !frequency.empty() || !time.empty(); }
	/// the material indices over each box, collective on the box communicators
	void write_masks();
	/// accumulates the frequency boxes at time t; writes the time boxes due at step
	void sample(double t, int step);
	/// writes the frequency boxes
	void write_dft();
};

CustomBoxes &custom_boxes();
} // namespace fdtd::monitors

#endif
