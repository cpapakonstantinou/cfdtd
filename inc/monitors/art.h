/**
 * @file monitors/art.h
 * @brief Absorbance, reflectance and transmittance of periodic structures from the ART planes
 *
 * @verbatim
      R = flux back through the reflection plane / incident,  T = flux through the transmission plane / incident
      A = 1 - R - T, and with [output] sigma_absorption the conductivity-based absorption
   @endverbatim
 * Written by root to reports/art_<n>.txt; n numbers the intermediate spectra (count_art).
 */
#ifndef FDTD_MONITORS_ART_H
#define FDTD_MONITORS_ART_H

namespace fdtd::monitors
{
/// collective, reductions over MPI_COMM_WORLD
void write_art_spectra(int step); ///< step: the time step the spectra are taken at, file names
/// this rank's conductivity-based absorption of material per frequency, sigma |E|^2 over its cells
/// into abs_sigma_loc
void integrate_absorption(int material);
} // namespace fdtd::monitors

#endif
