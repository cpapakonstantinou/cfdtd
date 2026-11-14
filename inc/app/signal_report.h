/** \file   inc/app/signal_report.h
 *  \brief  `fdtd signal`: the length of the source signal, the time steps a run needs and the
 *          spectrum it covers. No MPI, no fields.
 **/
#pragma once

namespace fdtd::app
{
/// runs the report on the arguments after the command, argv[0] being the command itself
int signal_report(int argc, char *argv[]);
} // namespace fdtd::app
