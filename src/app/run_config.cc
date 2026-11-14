/**
 * @file app/run_config.cc
 * @brief The run configuration object and its seal
 */
#include "app/run_config.h"
#include "app/run_state.h"
#include <stdexcept>

namespace fdtd::app
{
RunConfig the_run_config;
RunState the_run_state;

namespace
{
bool sealed = false;
} // namespace

RunConfig &edit_run_config()
{
	if (sealed)
		throw std::logic_error("run configuration changed after the configure phase");
	return the_run_config;
}

void seal_run_config() { sealed = true; }
void reopen_run_config() { sealed = false; }
} // namespace fdtd::app
