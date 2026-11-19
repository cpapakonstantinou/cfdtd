/**
 * @file grid/stretching.h
 * @brief The coordinate stretching of the absorbing slabs: 1 / kappa per axis, face and grading position
 *
 * The CPML stretches each coordinate u in its slabs, d/du -> (1 / kappa_u) d/du + psi_u. The kappa
 * part is a property of the coordinates: every derivative in a slab takes it, the field updates, the
 * hydrodynamic model, the 1D source grid, while the psi part belongs to the absorbing fields
 * boundary/cpml.h. The CPML fills these tables, Cpml::set_profiles.
 */
#ifndef FDTD_GRID_STRETCHING_H
#define FDTD_GRID_STRETCHING_H

#include "core/array.h"

namespace fdtd::grid
{
struct Stretching
{
	/// inv_kappa[axis][face][position]: position 0 the E samples in m, 1 the H samples (m + 1/2),
	/// m = 0, outer end .. N, interface
	array1<double> inv_kappa[3][2][2];
};
extern Stretching the_stretching;
inline Stretching &stretching() { return the_stretching; }
} // namespace fdtd::grid

#endif
