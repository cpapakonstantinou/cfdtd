/**
 * @file monitors/flux_boxes.h
 * @brief The DFT planes of the two flux boxes, scattering and total/absorption cross sections
 *
 * Each box has six faces; on each face the running DFT of the two tangential E and H components is
 * kept per frequency, for the field and, hydrodynamic runs, the harmonic field:
 * @verbatim
      box(b).face(n, side).re/im[field][set][c]     b: 0 scattering, 1 total
                                                    n: normal axis, side: 0 near, 1 far
                                                    field: 0 E, 1 H;  set: 0 field, 1 harmonic
                                                    c: component, only the two tangential to the face
      each array: [frequency][x][y][z] with size 0 along n
   @endverbatim
 * Read by the Q coefficients and ART-TFSF estimate (q_coefficients.cc) and the near-to-far faces
 * (near_to_far.cc).
 */
#ifndef FDTD_MONITORS_FLUX_BOXES_H
#define FDTD_MONITORS_FLUX_BOXES_H

#include "core/carray.h"

namespace fdtd::monitors
{
struct FluxFace
{
	ctetrad<double> re[2][2][3], im[2][2][3]; ///< [E/H][field/harmonic][component]
};

struct FluxBox
{
	FluxFace face[3][2]; ///< [normal axis][near/far]
};

struct FluxBoxes
{
	FluxBox box[2]; ///< scattering, total

	/// allocates the faces this rank has, q_coeff runs, harmonic sets with hydrodynamics
	void allocate();
	void release();
	/// accumulates this step: H at t_h, E at t_e
	void sample(double t_e, double t_h);
	/// divides the field set, not the harmonic set, by the source spectrum after the time loop
	void normalize();
};

FluxBoxes &flux_boxes();
} // namespace fdtd::monitors

#endif
