/**
 * @file grid/output_layout.h
 * @brief Where this rank's part of the physical domain, the grid without CPML, sits in the output files
 *
 * @verbatim
      start[a]   global index of this rank's first output cell
      global[a]  cells of the physical domain, total minus both CPML layers
      local[a]   this rank's inner cells, subdomain.inner_dim
      shift[a]   offset of the inner region in the rank's padded arrays, subdomain.inner_org
   @endverbatim
 * The binary writers (save_as_binary_2d / 3D) take these as int arrays.
 */
#ifndef FDTD_GRID_OUTPUT_LAYOUT_H
#define FDTD_GRID_OUTPUT_LAYOUT_H

namespace fdtd::grid
{
struct OutputLayout
{
	int start[3];
	int global[3];
	int local[3];
	int shift[3];
};

/// this rank's layout, valid once the domain is decomposed and the rank initialised
OutputLayout output_layout();
} // namespace fdtd::grid

#endif
