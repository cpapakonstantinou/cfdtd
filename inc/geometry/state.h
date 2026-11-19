/**
 * @file geometry/state.h
 * @brief State of the geometry module
 *
 * Moved out of the former shared namespace fdtd::state: read and written by this module, allocated or set by the
 * run set-up, app/setup.cc, and the configuration, applyConfig.cc. Defined in src/geometry/state.cc.
 */
#ifndef FDTD_GEOMETRY_STATE_H
#define FDTD_GEOMETRY_STATE_H

#include "core/array.h"
#include "core/constants.h"
#include "materials/material_field.h"
#include "materials/slots.h"
#include <complex>
#include <mpi.h>
#include <string>
#include <vector>

namespace fdtd::geometry
{
extern int j_offset_uni;                                    ///< Offset from the global y index to the 1D grid index used when filling background_y.
extern double corner_radius;                                  ///< Corner rounding of the rounded shapes in m;.
extern bool add_free_electron_layer;                        ///< Add a free-electron surface layer to the antenna.
extern double antenna_gap;                                            ///< Antenna gap in m.
} // namespace fdtd::geometry

#endif
