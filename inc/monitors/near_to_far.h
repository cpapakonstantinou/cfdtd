/**
 * @file monitors/near_to_far.h
 * @brief The scattering-box faces written for the Kirchhoff near-to-far transformation
 *
 * The normalised DFT of the tangential E and H on the six faces of the scattering box, the flux-box
 * faces of set 0, monitors/flux_boxes.h, one HDF5 file per face and wavelength:
 * @verbatim
      fields/<near|far>_<x|y|z>_<wl>nm.h5
          complex datasets E_a, E_b, H_a, H_b, a < b the axes tangential to the face
          attributes domain_size, face, wavelength
      fields/ini_near_to_far.dat, root, grid, box and centre for the transformation
   @endverbatim
 */
#ifndef FDTD_MONITORS_NEAR_TO_FAR_H
#define FDTD_MONITORS_NEAR_TO_FAR_H

namespace fdtd::monitors
{
/// writes the parameter file, root, and the faces, collective on the face communicators, a barrier
/// after each normal axis
void write_near_to_far();
} // namespace fdtd::monitors

#endif
