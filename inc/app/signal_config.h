/**
 * @file app/signal_config.h
 * @brief The run's source signal from the configuration, [source] signal, [signal.*], the band
 */
#ifndef FDTD_APP_SIGNAL_CONFIG_H
#define FDTD_APP_SIGNAL_CONFIG_H

namespace fdtd::app
{
/// builds the signal and sets it as the run's, sources::set_run_signal; no MPI. Called by the run
/// set-up, `fdtd check` and `fdtd signal` after the time step is set, the signal keeps dt.
void configure_signal();
} // namespace fdtd::app

#endif
