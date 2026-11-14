/**
 * @file monitors/inputs.h
 * @brief The monitors the configuration declares one by one: probe points and monitor boxes
 */
#ifndef FDTD_MONITORS_INPUTS_H
#define FDTD_MONITORS_INPUTS_H

#include <string>

namespace fdtd::monitors
{
struct ProbePoint
{
	std::string name;
	int x = 0, y = 0, z = 0;
};

/// one [monitor.<name>] section: a box of cells recorded in time or frequency
struct MonitorBox
{
	std::string name;
	std::string type = "frequency";      ///< frequency, DFT of E | time, E every `step`
	int center_x = 0, center_y = 0, center_z = 0; ///< [cells], inner domain, PML excluded
	int size_x = 0, size_y = 0, size_z = 0;       ///< [cells]
	int step = 1;                        ///< time monitors: record every N steps
};
} // namespace fdtd::monitors

#endif
