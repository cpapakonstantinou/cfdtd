/**
 * @file app/planes.h
 * @brief The global positions of the boxes and planes of this run, from the run configuration
 */
#ifndef FDTD_APP_PLANES_H
#define FDTD_APP_PLANES_H

#include "grid/lattice.h"
#include "grid/planes.h"

namespace fdtd::app
{
/// TF/SF box, ART planes, flux boxes and plot planes in global cells, CPML included
grid::GlobalPlanes global_planes(const grid::Lattice &lat);
} // namespace fdtd::app

#endif
