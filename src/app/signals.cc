/** \file   app/signals.cc
 *  \brief  Graceful handling of SIGINT and SIGTERM, inc/app/signals.h.
 */
#include "app/signals.h"

#include <csignal>

namespace fdtd::app
{
namespace
{
volatile std::sig_atomic_t received = 0; ///< the last signal received, 0 if none

extern void on_signal(int sig)
{
	received = sig;
	std::signal(sig, SIG_DFL); // a second signal ends the process at once
}
} // namespace

void install_signal_handlers()
{
	std::signal(SIGINT, on_signal);
	std::signal(SIGTERM, on_signal);
	std::signal(SIGUSR1, on_signal); // e.g. sbatch --signal=USR1@120: stop and dump before the time limit
}

int pending_signal() noexcept { return received; }

int agree_signal(MPI_Comm comm)
{
	int local = received;
	int any = 0;
	MPI_Allreduce(&local, &any, 1, MPI_INT, MPI_MAX, comm);
	return any;
}
} // namespace fdtd::app
