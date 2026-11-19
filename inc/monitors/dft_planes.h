/**
 * @file monitors/dft_planes.h
 * @brief The display DFT planes: frequency-domain field maps on the animation planes
 *
 * @verbatim
      plane[n]   the plane with normal n through the animation point, i/j/subdomain.k_animation:
                 re/im[E/H][field/harmonic][component], all three components
                 x: subdomain.plots_yz() && dft_plane_yz;  y: subdomain.plots_xz() && (dft_plane_xz || graphene_optim);  z: subdomain.plots_xy() && dft_plane_xy
      xz_1, xz_2 the backward / forward far-field xz planes, subdomain.j_animation_1 / _2: E only
      graphene   the xz plane one cell below subdomain.j_animation (graphene_optim): E only
      each array: [frequency][x][y][z] with size 0 along the normal
   @endverbatim
 */
#ifndef FDTD_MONITORS_DFT_PLANES_H
#define FDTD_MONITORS_DFT_PLANES_H

#include "core/array.h"

namespace fdtd::monitors
{
struct DftPlane
{
	array4<double> re[2][2][3]; ///< [E/H][field/harmonic][component]
	array4<double> im[2][2][3]; ///< [E/H][field/harmonic][component]
};

struct DftPlanes
{
	DftPlane plane[3];               ///< normal x (yz), y (xz), z (xy)
	DftPlane xz_1; ///< E only (re/im[0][0][c])
	DftPlane xz_2; ///< E only (re/im[0][0][c])
	DftPlane graphene; ///< E only (re/im[0][0][c])

	/// allocates the planes this rank has, not on ranks that are all CPML
	void allocate();
	void release();
	/// accumulates this step at time t, E and H alike
	void sample(double t);
	/// divides by the source spectrum after the time loop, the field set; not the harmonic planes
	void normalize();
	/// writes DFT/<E|H>[_nonlinear|_graphene|_forward|_backward]_<wl>nm_<XY|XZ|YZ>.h5: complex datasets
	/// x, y, z, attributes domain_size, plane, wavelength, collective on the plot communicators
	void write_maps();
};

DftPlanes &dft_planes();
} // namespace fdtd::monitors

#endif
