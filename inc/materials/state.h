/**
 * @file materials/state.h
 * @brief State of the materials module
 *
 * Moved out of the former shared namespace fdtd::state: read and written by this module, allocated or set by the
 * run set-up, app/setup.cc, and the configuration, applyConfig.cc. Defined in src/materials/state.cc.
 */
#ifndef FDTD_MATERIALS_STATE_H
#define FDTD_MATERIALS_STATE_H

#include "core/carray.h"
#include "core/constants.h"
#include "materials/material_field.h"
#include "materials/material_map.h"
#include "materials/codes.h"
#include "materials/dispersive_state.h"
#include <complex>
#include <mpi.h>
#include <string>
#include <vector>

namespace fdtd::materials
{
extern cvector<int> background_y; ///< background material code per global y cell, the layered medium the 1D plane-wave grid and the CPML faces see
extern bool laplacian;                                      ///< Nonlocal term uses the Laplacian of P.
extern bool gradient_divergence;                            ///< Nonlocal term uses grad(div P).
extern ctensor<double> magnetization_x;                      ///< Static magnetisation x component read from file, per cell, magnetic Lorentz.
extern ctensor<double> magnetization_y;                      ///< Static magnetisation y component read from file, per cell, magnetic Lorentz.
extern ctensor<double> magnetization_z;                      ///< Static magnetisation z component read from file, per cell, magnetic Lorentz.
} // namespace fdtd::materials

#endif
