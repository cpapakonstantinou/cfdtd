/**
 * @file materials/thermal/state.cc
 * @brief Definitions of the thermal module state, thermal/state.h
 */
#include "materials/thermal/state.h"

namespace fdtd::materials::thermal
{
double pump_probe_delay;
int nr_max_time_steps_pump;
std::string str_outfilepathname_t;
std::string str_infilepathname_t;
std::string str_fileextension;
int num_rk;
fdtd::array3<double> electron_temp;
fdtd::array3<double> electron_temp_tf;
fdtd::array3<double> lattice_temp;
fdtd::array3<double> temp_save;
} // namespace fdtd::materials::thermal
