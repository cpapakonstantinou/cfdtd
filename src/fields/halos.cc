/**
 * @file fields/halos.cc
 * @brief Halo exchange of the fields object, fields/halos.h, over the grid's exchange primitive
 */
#include "fields/halos.h"
#include "fields/fields.h"
#include "grid/state.h"
#include "grid/subdomain.h"

namespace fdtd::fields
{
using grid::exchange;
using grid::Field;
using grid::subdomain;
using grid::xy_interface_data_type;
using grid::xz_interface_data_type;
using grid::yz_interface_data_type;

void exchange_halos(Field field, int set)
{
	auto &F = fdtd::fields::fields();
	auto &c = field == Field::E ? F.E[set] : F.H[set];
	const auto &S = subdomain;
	exchange(field, 0, c[1], c[2], S.neighbours_x[0], S.neighbours_x[1], S.cells_x, yz_interface_data_type);
	exchange(field, 1, c[0], c[2], S.neighbours_y[0], S.neighbours_y[1], S.cells_y, xz_interface_data_type);
	exchange(field, 2, c[0], c[1], S.neighbours_z[0], S.neighbours_z[1], S.cells_z, xy_interface_data_type);
}

void exchange_longitudinal_e()
{
	auto &E = fdtd::fields::fields().E[0];
	const auto &S = subdomain;
	exchange(Field::E, 0, E[0], E[0], S.neighbours_nl_x[0], S.neighbours_nl_x[1], S.cells_x, yz_interface_data_type);
	exchange(Field::E, 1, E[1], E[1], S.neighbours_nl_y[0], S.neighbours_nl_y[1], S.cells_y, xz_interface_data_type);
	exchange(Field::E, 2, E[2], E[2], S.neighbours_nl_z[0], S.neighbours_nl_z[1], S.cells_z, xy_interface_data_type);
}
} // namespace fdtd::fields
