/**
 * @file monitors/volume_dft.h
 * @brief The volume DFT: running DFT of E, and H, over a box of each rank's inner cells
 *
 * @verbatim
      field[E/H][c]   the running DFT, normalised in place after the time loop, VolumeDft::normalize;
                      E for the 3D outputs, the absorption integral or the multipoles; H with the 3D outputs
      pnl[c]          DFT of the nonlinear polarisation over the harmonic band, [m - harmonic_start],
                      exact multipoles with SHG/THG and more than two frequencies
      each array:     [frequency][x][y][z] over the cells of the box
   @endverbatim
 *
 * The box is the whole inner domain, with its ghost plane, when an output needs the whole volume: the
 * 3D and plane field saves or the absorption integral. For the multipoles alone it is the bounding box
 * of the cells the moments integrate on this rank, monitors/exact_multipoles.h, multipole_box; a rank
 * without such cells holds no arrays.
 */
#ifndef FDTD_MONITORS_VOLUME_DFT_H
#define FDTD_MONITORS_VOLUME_DFT_H

#include "core/array.h"
#include <complex>

namespace fdtd::monitors
{
/// a box of local inner cells: n cells per axis from lo
struct CellBox
{
	int lo[3] = {0, 0, 0};
	int n[3] = {0, 0, 0};

	bool empty() const { return n[0] <= 0 || n[1] <= 0 || n[2] <= 0; }
	bool contains(int i, int j, int k) const
	{
		return i >= lo[0] && i < lo[0] + n[0] && j >= lo[1] && j < lo[1] + n[1] && k >= lo[2] && k < lo[2] + n[2];
	}
};

struct VolumeDft
{
	array4<std::complex<double>> field[2][3];
	array4<std::complex<double>> pnl[3];
	CellBox box;                 ///< the cells the arrays cover
	int harmonic_start = 0;      ///< the first frequency of pnl

	/// after the mesh: the box this run needs and the arrays over it, none on ranks that are all CPML;
	/// false when the memory runs out
	bool allocate();
	/// releases the arrays
	void release();
	/// one step of the running DFT of E at time t
	void sample_e(double t);
	/// one step of the running DFT of H at time t
	void sample_h(double t);
	/// field[f] divided by S(w_m), f: 0 E, 1 H
	void normalize(int f);
	/// pnl divided by S(w_m), SHG/THG with more than two frequencies
	void normalize_pnl();

	/// field f, component c, frequency m at local inner cell (i, j, k) of the box
	std::complex<double> &at(int f, int c, int m, int i, int j, int k) { return field[f][c][m][i - box.lo[0]][j - box.lo[1]][k - box.lo[2]]; }
	const std::complex<double> &at(int f, int c, int m, int i, int j, int k) const
	{
		return field[f][c][m][i - box.lo[0]][j - box.lo[1]][k - box.lo[2]];
	}
	/// the nonlinear polarisation, component c, frequency m of the harmonic band, at local inner cell (i, j, k) of the box
	std::complex<double> &pnl_at(int c, int m, int i, int j, int k) { return pnl[c][m - harmonic_start][i - box.lo[0]][j - box.lo[1]][k - box.lo[2]]; }
};

VolumeDft &volume_dft();
} // namespace fdtd::monitors

#endif
