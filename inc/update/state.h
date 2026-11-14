/**
 * @file update/state.h
 * @brief State of the update module
 *
 * Moved out of the former shared namespace fdtd::state: read and written by this module, allocated or set by the
 * run set-up, app/setup.cc, and the configuration, applyConfig.cc. Defined in src/update/state.cc.
 */
#ifndef FDTD_UPDATE_STATE_H
#define FDTD_UPDATE_STATE_H

#include "core/carray.h"
#include "core/constants.h"
#include "materials/material_field.h"
#include "materials/codes.h"
#include <complex>
#include <mpi.h>
#include <string>
#include <vector>

namespace fdtd::update
{
} // namespace fdtd::update

#endif
