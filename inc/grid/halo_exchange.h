/** \file   grid/halo_exchange.h
 *  \brief  Halo exchange of the fields between neighbouring ranks, src/grid/halo_exchange.cc.
 */
#ifndef FDTD_GRID_HALO_EXCHANGE_H
#define FDTD_GRID_HALO_EXCHANGE_H

#include "mpi.h"
#include "core/views.h"

namespace fdtd::grid
{
/// which way a field's halo is filled: E from the lower neighbour, H from the upper one
enum class Field { E, H };

/**
 * @brief Exchanges the halo planes normal to `axis` (0 x, 1 y, 2 z) of the two field components
 *        f1, f2 with the neighbours less / more; n cells along the axis; `plane` the MPI type of one
 *        plane. Handles the periodic boundaries along x and z.
 */
void exchange(Field field, int axis, view3<double> f1, view3<double> f2, int less, int more, int n, MPI_Datatype plane);
// the exchanges of the fields object are fields::exchange_halos / exchange_longitudinal_e, fields/halos.h
/// the margin exchange of the hydrodynamic model and the LWA multipoles, two or three planes per
/// component and face; irregular, to be folded with the hydrodynamic review: three components / one
void exchange_margins(view3<double> x, view3<double> y, view3<double> z);
void exchange_margin(view3<double> x);
} // namespace fdtd::grid

#endif // FDTD_GRID_HALO_EXCHANGE_H
