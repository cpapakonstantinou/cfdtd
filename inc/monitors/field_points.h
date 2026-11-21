/**
 * @file monitors/field_points.h
 * @brief Ex and Ez, amplitude, phase, ratio, phase difference, per wavelength at single points of
 *        the normalised DFT planes, for the antenna shapes
 *
 * @verbatim
      far_field_forward.txt         periodic ART runs: the transmission plane;  SHG/THG cylinder: the forward xz plane
      far_field_backward.txt        periodic ART runs: the reflection plane;    SHG/THG cylinder: the backward xz plane
      near_field.txt                the xz display plane, else the yz plane, at the animation point
      each at, subdomain.i_animation, subdomain.k_animation, of the plane; columns LAMBDA[nm], |Ex|, |Ez|, |Ex|/|Ez|,
      phase Ex, phase Ez, phase difference, degrees, wrapped to [-180, 180]
   @endverbatim
 */
#ifndef FDTD_MONITORS_FIELD_POINTS_H
#define FDTD_MONITORS_FIELD_POINTS_H

namespace fdtd::monitors
{
/// writes the point files this rank's planes and the shape call for, rank-local
void write_field_points();
} // namespace fdtd::monitors

#endif
