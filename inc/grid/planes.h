/**
 * @file grid/planes.h
 * @brief The boxes and planes the decomposition maps onto the ranks, in global cells, CPML included
 *
 * The grid module's second input next to the Lattice: the application decides where the TF/SF box,
 * the ART planes, the flux boxes and the plot planes are, app::global_planes; set_lattice_data finds
 * which ranks they cross and their local indices.
 */
#ifndef FDTD_GRID_PLANES_H
#define FDTD_GRID_PLANES_H

namespace fdtd::grid
{
struct GlobalPlanes
{
	int tfsf_nea[3] = {}; ///< TF/SF box faces along x, y, z
	int tfsf_far[3] = {}; ///< TF/SF box faces along x, y, z
	int art_nea_y = 0; ///< ART reflection / transmission planes
	int art_far_y = 0; ///< ART reflection / transmission planes
	int qsca_nea[3] = {}; ///< scattering flux box faces
	int qsca_far[3] = {}; ///< scattering flux box faces
	int qtot_nea[3] = {}; ///< total flux box faces
	int qtot_far[3] = {}; ///< total flux box faces
	int plot_x = 0; ///< plot planes (yz, xz, xy)
	int plot_y = 0; ///< plot planes (yz, xz, xy)
	int plot_z = 0; ///< plot planes (yz, xz, xy)
	int plot_y1 = 0; ///< the two extra xz plot planes
	int plot_y2 = 0; ///< the two extra xz plot planes
	bool plot_y_not_on_rank_start = false;  ///< fail if the xz plot plane is the first cell of a rank (graphene_optim)
};
} // namespace fdtd::grid

#endif
