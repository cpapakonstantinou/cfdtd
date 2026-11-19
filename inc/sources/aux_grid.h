/**
 * @file sources/aux_grid.h
 * @brief The 1D auxiliary grid of the TF/SF plane-wave source
 *
 * The plane wave is propagated along y on a 1D grid of the same layered media, first medium,
 * second medium, and sampled on the TF/SF faces of the 3D grid. Two polarisations run side by side:
 * @verbatim
      pol 0 (x): Ex, Hz      curl Ex = dHz/dy        Hz += dt/mu0 dEx/dy
      pol 1 (z): Ez, Hx      curl Ez = -dHx/dy       Hx -= dt/mu0 dEz/dy
   @endverbatim
 * Each polarisation has its own ADE state, Drude/Lorentz pole, critical points, and CPML psi at the
 * two ends. Line 0 carries the source; line 1 is a second line in the first medium only, kept for periodic
 * x and z, where it gives the incident field on the near y face.
 */
#ifndef FDTD_SOURCES_AUX_GRID_H
#define FDTD_SOURCES_AUX_GRID_H

#include "core/array.h"
#include "materials/slots.h"
#include "materials/state.h"

namespace fdtd::sources
{
/// one 1D line along y; index p is the polarisation (0 x, 1 z), side 0 near / 1 far
struct AuxLine
{
	array1<double> e[2];       ///< Ex, Ez
	array1<double> e_prev[2];  ///< E one step back
	array1<double> h[2];       ///< Hz, Hx, the H of each polarisation
	array1<double> p[2];       ///< Drude / Lorentz polarisation
	array1<double> p_prev[2];  ///< and one step back
	array1<double> pc[2][materials::max_poles];      ///< critical-point polarisations per pole
	array1<double> pc_prev[2][materials::max_poles]; ///< and one step back
	array1<double> psi_e[2][2]; ///< CPML psi of E (d/dy), [p][side], num_cpml_cells_y + 1 cells
	array1<double> psi_h[2][2]; ///< CPML psi of H (d/dy), [p][side]

	/// allocates and zeroes every array: fields on n_cells + 1 points, psi on n_cpml + 1
	void allocate(int n_cells, int n_cpml);
	void release();
};

/// the auxiliary grid: the lines and the material of each 1D cell
struct AuxGrid
{
	AuxLine line[2];     ///< 0 the source line, 1 the periodic reference line
	array1<int> &mat = fdtd::materials::background_y; ///< material slot per cell of line 0, the background along y

	/// E update of line 0 and, with periodic x and z, line 1, the TF/SF source step before the 3D E scan
	void update_e();
	/// H update of the same lines
	void update_h();
	void release();
};

/// the auxiliary grid of this rank
AuxGrid &aux_grid();
} // namespace fdtd::sources

#endif
