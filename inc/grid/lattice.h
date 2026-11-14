/**
 * @file grid/lattice.h
 * @brief The global lattice the grid module decomposes: cells, CPML, periodicity, cell size.
 *
 * The grid module's input. The application sets it once from the run configuration, after the
 * periodic plane-wave adjustment of the CPML, app::setup_domain, and before the decomposition; the
 * grid module reads it instead of the run configuration, docs/dev/Roadmap.md §2.5.
 */
#ifndef FDTD_GRID_LATTICE_H
#define FDTD_GRID_LATTICE_H

namespace fdtd::grid
{
struct Lattice
{
	int cells[3] = {};      ///< global cells along x, y, z, CPML included
	int cpml[3] = {};       ///< CPML cells on each face along x, y, z
	int tfsf_inset[3] = {}; ///< TF/SF box inset from the inner domain along x, y, z in cells
	bool periodic_x = false, periodic_z = false; ///< periodic boundaries along x, z
	double d[3] = {};       ///< cell size along x, y, z in m
	double inv_d[3] = {};   ///< 1 / d in 1/m; set by set_lattice
	int inner(int a) const { return cells[a] - 2 * cpml[a]; } ///< global cells along a without the CPML
};

/// the lattice, valid once the application has set it
const Lattice &lattice();
/// set by the application before the decomposition, computes inv_d
void set_lattice(const Lattice &l);
} // namespace fdtd::grid

#endif
