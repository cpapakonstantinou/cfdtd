/**
 * @file monitors/exact_multipoles.h
 * @brief The exact multipole decomposition, Alaee, Rockstuhl & Fernandez-Corbaton, Opt. Commun.
 *        407, 17, 2018, of the induced current inside the antenna, from the normalised volume DFT
 *
 * @verbatim
      J_w(r) = -i w eps0 [eps_r(r, w) - eps_host(w)] E_w(r)          contrast current, all frequencies
      J_tot  = -i w [P_eq + P_NL,w]                                  harmonic band, SHG/THG: plus the
                                                                     DFT of the nonlinear polarisation
      moments p, m, Qe, Qm with the spherical-Bessel weights of Table 2, origin at the antenna centre
   @endverbatim
 * Output, root, reports/: multipoles_exact.txt and, with SHG/THG, multipoles_exact_full_current.txt,
 * each with its _cross_sections.txt; multipoles_lwa*.txt likewise. Details in the source.
 */
#ifndef FDTD_MONITORS_EXACT_MULTIPOLES_H
#define FDTD_MONITORS_EXACT_MULTIPOLES_H

#include "monitors/volume_dft.h"

namespace fdtd::monitors
{
/// the bounding box of this rank's inner cells the moments integrate: a Yee component in the
/// multipole volume or, with SHG/THG and a harmonic band, of a nonlinear material; empty if none
CellBox multipole_box();
/// one step of the DFT of P_NL over the harmonic band, pnl of the volume DFT, at time t
void sample_nonlinear_polarisation(double t);
/// the moments and cross sections; collective, ranks with only CPML contribute zeros; returns at
/// once unless [multipoles] exact is on
void write_exact_multipoles();
} // namespace fdtd::monitors

#endif
