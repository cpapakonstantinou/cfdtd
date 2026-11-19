/**
 * @file monitors/plane_dft.h
 * @brief The running DFT on one plane of the inner region and its normalisation, display planes,
 *        flux faces, ART planes
 *
 * @verbatim
      sample:     re[c][m][p] += f_c(p) cos(h w_m t),  im[c][m][p] += f_c(p) sin(h w_m t)
                  h = 1 for the field set, 2 for the harmonic set, the harmonic fields when they are
                  kept separately, else the fields
   @endverbatim
 *
 * @verbatim
      phase_normalization:  (re + i im) <- (re + i im) / S(w_m)     S: complex DFT of the source signal
      otherwise:            re <- re / |S(w_m)|,  im <- im / |S(w_m)|
      over [frequency][x][y][z] with index 0 along the normal and 0..inner_dim inclusive across it
   @endverbatim
 * Components without an array, the normal component of a flux face, E_y of the ART planes, are skipped.
 */
#ifndef FDTD_MONITORS_PLANE_DFT_H
#define FDTD_MONITORS_PLANE_DFT_H

#include "core/array.h"

namespace fdtd::monitors
{
/// one step of the DFT of E, field 0, or H, field 1, of set `set` at time t on the plane with normal
/// axis `normal` at `pad` cells from the inner origin; re/im[component], components without an array skipped
void sample_plane(int normal, int pad, int field, int set, double t, const array4<double> (&re)[3], const array4<double> (&im)[3]);
/// normalises the plane with normal axis `normal` in place; re/im[component]
void normalize_plane(int normal, const array4<double> (&re)[3], const array4<double> (&im)[3]);
} // namespace fdtd::monitors

#endif
