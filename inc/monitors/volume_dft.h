/**
 * @file monitors/volume_dft.h
 * @brief The volume DFT: running DFT of E, and H, over the whole inner domain
 *
 * @verbatim
      re/im[E/H][c]      running DFT; E with save_output, absorption (ABS_SIGMA), DFG optimisation or
                         the exact multipoles; H with save_output
      field[E/H][c]      the normalised complex DFT, normalizeDFT.cc, written as the 3D / plot-plane
                         output; allocated with E
      re/im_pnl[c]       DFT of the nonlinear polarisation over the harmonic band, exact multipoles
                         with SHG/THG and more than two frequencies
      each array: [frequency][x][y][z] over the inner cells
   @endverbatim
 */
#ifndef FDTD_MONITORS_VOLUME_DFT_H
#define FDTD_MONITORS_VOLUME_DFT_H

#include "core/carray.h"
#include <complex>

namespace fdtd::monitors
{
struct VolumeDft
{
	ctetrad<double> re[2][3], im[2][3];
	ctetrad<std::complex<double>> field[2][3];
	ctetrad<double> re_pnl[3], im_pnl[3];

	/// allocates what the run needs, not on ranks that are all CPML; arrays already there are kept
	void allocate();
	/// releases the arrays
	void release();
	/// one step of the running DFT of E at time t (re/im[0])
	void sample_e(double t);
	/// one step of the running DFT of H at time t (re/im[1])
	void sample_h(double t);
	/// field[f] = (re + i im)[f] / S(w_m), and re/im[f] replaced by it (f: 0 E, 1 H)
	void normalize(int f);
	/// re/im_pnl over the harmonic band divided by S(w_m), SHG/THG with more than two frequencies
	void normalize_pnl();
};

VolumeDft &volume_dft();
} // namespace fdtd::monitors

#endif
