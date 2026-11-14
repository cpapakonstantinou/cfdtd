/**
 * @file app/state_file.h
 * @brief The state dump of an interrupted run and its restart
 *
 * When a signal, or the wall-clock limit, stops the time loop, every rank writes
 * <output>/state/rank_<r>.bin: the step, the fields, the CPML and auxiliary-grid state, the dispersive
 * state and the monitors' accumulators. `fdtd simulate --restart` sets the run up as usual, loads the
 * files and continues the time loop at the next step; the result is the one of an uninterrupted run.
 * The restart must use the same configuration and rank layout, checked.
 */
#ifndef FDTD_APP_STATE_FILE_H
#define FDTD_APP_STATE_FILE_H

#include <string>

namespace fdtd::app
{
/// why this run can not be dumped and restarted, empty when it can
std::string state_unsupported();
/// writes this rank's state file after `step` completed steps; false, and `why`, on failure
bool save_state(const std::string &dir, int step, std::string &why);
/// loads this rank's state file into the allocated run; the completed steps in `step`
bool load_state(const std::string &dir, int &step, std::string &why);
} // namespace fdtd::app

#endif
