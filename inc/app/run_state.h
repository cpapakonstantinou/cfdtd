/**
 * @file app/run_state.h
 * @brief The state of the run that changes while it runs: time and step
 *
 * Owned by the driver, app; the time loop writes it, the modules read it.
 * The run configuration, fixed for the run, is app::RunConfig.
 */
#ifndef FDTD_APP_RUN_STATE_H
#define FDTD_APP_RUN_STATE_H

#include <cstdio>
#include <ctime>

namespace fdtd::app
{
/// time and step of the run
struct RunState
{
	int nr_max_time_steps{};                                ///< Number of time steps the run will take.
	double t{};                                             ///< Current simulation time in s. A global named t; avoid shadowing it.
	double t_1{};                                           ///< Time of the H fields this step, t + DT/2 in s.
	int timestep{};                                         ///< Current time step index.
	clock_t prgstart{};                                     ///< Clock at program start, for log timestamps.
	FILE * ex_file = nullptr;                               ///< Text file receiving the Ex probe trace (Ex.txt).
	FILE * ez_file = nullptr;                               ///< Text file receiving the Ez probe trace (Ez.txt).
	char cstring[100]{};                                    ///< Scratch buffer for formatting step numbers.
};

extern RunState the_run_state;
/// the run state
inline RunState &run_state() { return the_run_state; }
} // namespace fdtd::app

#endif
