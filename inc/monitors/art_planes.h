/**
 * @file monitors/art_planes.h
 * @brief The reflection and transmission planes of periodic runs, ART: absorbance, reflectance,
 *        transmittance
 *
 * Two xz planes, each with the running DFT of Ex, Ez, Hx, Hz:
 * @verbatim
      plane[0]  reflection, near y, subdomain.j_art_1, on ranks with art_nea_y
      plane[1]  transmission, far y, subdomain.j_art_2, on ranks with art_far_y
      re/im[E/H][0][component], components x and z, [1] stays empty; [frequency][x][0][z]
   @endverbatim
 * Before 2026-11-07 both planes accumulated into one set of arrays, so a rank holding both, every run
 * with one rank along y, integrated R and T over their sum (T = -R).
 */
#ifndef FDTD_MONITORS_ART_PLANES_H
#define FDTD_MONITORS_ART_PLANES_H

#include "monitors/dft_planes.h"

namespace fdtd::monitors
{
struct ArtPlanes
{
	DftPlane plane[2]; ///< reflection, near, transmission, far

	bool has(int side) const;
	void allocate();
	void release();
	/// accumulates this step: E at t_e, H at t_h
	void sample(double t_e, double t_h);
	/// divides by the DFT normalisation, after the time loop
	void normalize();
	/**
	 * @brief Poynting flux through plane side per frequency (out[m]), +y outward: Re(Ez Hx* - Ex Hz*)
	 * @return the plane's area on this rank
	 */
	double flux(int side, double out[]) const;
};

ArtPlanes &art_planes();
} // namespace fdtd::monitors

#endif
