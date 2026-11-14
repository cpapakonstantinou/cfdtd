/**
 * @file monitors/q_coefficients.h
 * @brief Scattering, absorption and extinction efficiencies from the flux boxes, and the TF/SF
 *        estimate of reflectance and transmittance
 *
 * @verbatim
      Q_sca = (outward flux through the scattering box) / (incident intensity x area)
      Q_abs = (inward flux through the total-field box) / (incident intensity x area)
      area: the shape's geometric cross section, else the TF/SF box face
   @endverbatim
 * Written by root to ART/Qcoeff_<prec>_<pol>.txt, and the harmonic NL_Efficiency files, and
 * reports/ART-TFSF_<prec>_<pol>_<n>.txt.
 */
#ifndef FDTD_MONITORS_Q_COEFFICIENTS_H
#define FDTD_MONITORS_Q_COEFFICIENTS_H

namespace fdtd::monitors
{
/// the efficiencies; collective
void write_q_coefficients();
/// R and T approximated by single faces of the boxes, non-periodic runs; collective
void write_art_tfsf();
} // namespace fdtd::monitors

#endif
