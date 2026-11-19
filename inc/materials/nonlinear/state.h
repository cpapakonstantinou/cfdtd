/**
 * @file materials/nonlinear/state.h
 * @brief State of the nonlinear module
 *
 * Moved out of the former shared namespace fdtd::state: read and written by this module, allocated or set by the
 * run set-up, app/setup.cc, and the configuration, applyConfig.cc. Defined in src/materials/nonlinear/state.cc.
 */
#ifndef FDTD_MATERIALS_NONLINEAR_STATE_H
#define FDTD_MATERIALS_NONLINEAR_STATE_H

#include "core/array.h"
#include "core/constants.h"
#include "materials/material_field.h"
#include "materials/slots.h"
#include <complex>
#include <mpi.h>
#include <string>
#include <vector>

namespace fdtd::materials::nonlinear
{
extern double chi_2[materials::max_slots][18]; ///< chi(2) = 2 d, 18 Voigt entries, per slot; zero outside SHG runs
extern double chi_3[materials::max_slots][21]; ///< chi(3), 21 entries, per slot; zero outside THG runs
extern bool iso_chi_2[materials::max_slots];   ///< chi(2) of the slot has diagonal entries only
extern bool iso_chi_3[materials::max_slots];   ///< chi(3) of the slot has diagonal entries only
extern bool nl[materials::max_slots];          ///< the slot is nonlinear in this run: a tensor entry is not zero
} // namespace fdtd::materials::nonlinear

#endif
