/**
 * @file materials/thermal/state.h
 * @brief State of the thermal module
 *
 * Moved out of the former shared namespace fdtd::state: read and written by this module, allocated or set by the
 * run set-up, app/setup.cc, and the configuration, applyConfig.cc. Defined in src/materials/thermal/state.cc.
 */
#ifndef FDTD_MATERIALS_THERMAL_STATE_H
#define FDTD_MATERIALS_THERMAL_STATE_H

#include "core/carray.h"
#include "core/constants.h"
#include "materials/material_field.h"
#include "materials/codes.h"
#include <complex>
#include <mpi.h>
#include <string>
#include <vector>

namespace fdtd::materials::thermal
{
extern double pump_probe_delay;                                 ///< Delay between pump and probe pulses in s.
extern int nr_max_time_steps_pump;                          ///< Number of steps of the pump run.
extern std::string str_outfilepathname_t;                        ///< Path prefix of the temperature output files.
extern std::string str_infilepathname_t;                         ///< Path prefix of the temperature input file.
extern std::string str_fileextension;                            ///< File extension of the temperature files.
extern int num_rk;                                          ///< Number of Runge-Kutta stages of the temperature update.
extern ctensor<double> electron_temp;                         ///< Electron temperature per cell in K.
extern ctensor<double> electron_temp_tf;                      ///< Electron temperature written to / read from the temperature files in K.
extern ctensor<double> lattice_temp;                          ///< Lattice temperature per cell in K.
extern ctensor<double> temp_save;                             ///< Temperature history of the pump run [field][step][y], replayed in the probe run.
} // namespace fdtd::materials::thermal

#endif
