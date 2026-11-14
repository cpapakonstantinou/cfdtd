/**
 * @file grid/lattice.cc
 * @brief The global lattice, grid/lattice.h, and this rank's output layout, grid/output_layout.h
 */
#include "grid/output_layout.h"
#include "grid/lattice.h"
#include "grid/stretching.h"
#include "grid/subdomain.h"

namespace fdtd::grid
{
OutputLayout output_layout()
{
	OutputLayout L{};
	for (int a = 0; a < 3; ++a)
		L.start[a] = int(subdomain.output_start[a]); // long -> int: the binary writers take int offsets
	L.global[0] = lattice().cells[0] - 2 * lattice().cpml[0];
	L.global[1] = lattice().cells[1] - 2 * lattice().cpml[1];
	L.global[2] = lattice().cells[2] - 2 * lattice().cpml[2];
	L.local[0] = subdomain.inner_dim_x;
	L.local[1] = subdomain.inner_dim_y;
	L.local[2] = subdomain.inner_dim_z;
	L.shift[0] = subdomain.inner_org_x;
	L.shift[1] = subdomain.inner_org_y;
	L.shift[2] = subdomain.inner_org_z;
	return L;
}

namespace
{
Lattice the_lattice; ///< set once by the application
} // namespace
const Lattice &lattice() { return the_lattice; }
void set_lattice(const Lattice &l)
{
	the_lattice = l;
	for (int a = 0; a < 3; ++a)
		the_lattice.inv_d[a] = double(1.) / l.d[a];
}
} // namespace fdtd::grid

fdtd::grid::Stretching fdtd::grid::the_stretching;
