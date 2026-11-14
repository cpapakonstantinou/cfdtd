/** \file   app/signals.h
 *  \brief  Graceful handling of SIGINT and SIGTERM.
 *
 *  The handler only records the signal, a std::sig_atomic_t, the one thing a handler may safely
 *  write; see https://en.cppreference.com/cpp/utility/program/signal, and restores the default
 *  action, so a second signal ends the process at once. The time loop asks agree_signal() once
 *  per step: a signal may reach only some ranks, so the ranks agree on it, MPI_Allreduce, MAX
 *  and all leave the loop at the same step, with no rank left waiting in a halo exchange. The
 *  state machine then finishes the sweep, writes the outputs and exits with 128 + signal.
 */
#ifndef FDTD_APP_SIGNALS_H
#define FDTD_APP_SIGNALS_H

#include "mpi.h"

namespace fdtd::app
{
/// installs the handler for SIGINT and SIGTERM
void install_signal_handlers();

/// the signal this rank received, or 0
int pending_signal() noexcept;

/// the largest signal number received by any rank of comm, or 0, collective
int agree_signal(MPI_Comm comm);
} // namespace fdtd::app

#endif // FDTD_APP_SIGNALS_H
