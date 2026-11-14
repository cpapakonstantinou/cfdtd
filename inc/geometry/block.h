/**
 * @file geometry/block.h
 * @brief A box of one material and its copies, [block.<name>] in the configuration
 */
#ifndef FDTD_GEOMETRY_BLOCK_H
#define FDTD_GEOMETRY_BLOCK_H

#include <string>
#include <vector>

namespace fdtd::geometry
{
/** one [block.<name>] section: a box of one material, drawn before the [geometry] shape, in file order
 *  a later block overwrites an earlier one. Bounds are from the structure centre; the default
 *  +-1 m reaches past the domain on that side. `count` copies at `pitch_*` make a grating. */
struct Block
{
	std::string name;
	int material = 0;                       ///< material code
	double x_min = -1, x_max = 1;           ///< [m]
	double y_min = -1, y_max = 1;           ///< [m]
	double z_min = -1, z_max = 1;           ///< [m]
	int count = 1;                          ///< copies
	double pitch_x = 0, pitch_y = 0, pitch_z = 0; ///< [m] between copies
};

/// draw the blocks in order, a later block overwrites an earlier one, about the scene's centre cell
void create_blocks(const std::vector<Block> &blocks);
} // namespace fdtd::geometry

#endif
