/**
 * @file materials/state.h
 * @brief State of the materials module
 *
 * Moved out of the former shared namespace fdtd::state: read and written by this module, allocated or set by the
 * run set-up, app/setup.cc, and the configuration, applyConfig.cc. Defined in src/materials/state.cc.
 */
#ifndef FDTD_MATERIALS_STATE_H
#define FDTD_MATERIALS_STATE_H

#include "core/array.h"
#include "core/constants.h"
#include "materials/material_field.h"
#include "materials/material_map.h"
#include "materials/slots.h"
#include "materials/dispersive_state.h"
#include <complex>
#include <mpi.h>
#include <string>
#include <vector>

namespace fdtd::materials
{
extern array1<int> background_y; ///< background material slot per global y cell, the layered medium the 1D plane-wave grid and the CPML faces see
extern bool laplacian;                                      ///< Nonlocal term uses the Laplacian of P.
extern bool gradient_divergence;                            ///< Nonlocal term uses grad(div P).
extern array3<double> magnetization_x;                      ///< Static magnetisation x component read from file, per cell, magnetic Lorentz.
extern array3<double> magnetization_y;                      ///< Static magnetisation y component read from file, per cell, magnetic Lorentz.
extern array3<double> magnetization_z;                      ///< Static magnetisation z component read from file, per cell, magnetic Lorentz.
} // namespace fdtd::materials

#endif
