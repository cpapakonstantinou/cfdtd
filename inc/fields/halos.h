/**
 * @file fields/halos.h
 * @brief Halo exchange of the E and H fields with the neighbour ranks
 */
#ifndef FDTD_FIELDS_HALOS_H
#define FDTD_FIELDS_HALOS_H

#include "grid/halo_exchange.h"

namespace fdtd::fields
{
/// the halos of one set of the fields object, 0 linear, 1 harmonic, of E or H along x, y and z: on
/// each axis the two transverse components, with the rank's neighbours
void exchange_halos(grid::Field field, int set);
/// each E component along its own axis, the longitudinal halo the nonlinear E update reads, with
/// the neighbours of the nonlinear layout
void exchange_longitudinal_e();
} // namespace fdtd::fields

#endif
