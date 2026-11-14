/**
 * @file sources/mode_source.h
 * @brief A waveguide mode from MPB, injected one way on a plane normal to x, y or z
 *
 * `[source] type = mode`. The mode's transverse E and H are read from MPB's field files of one band
 * a 2D cross-section, interpolated onto the Yee positions of the injection plane, and added as a
 * total-field / scattered-field correction on that plane: the total field is downstream, +axis, only
 * the scattered field upstream. With propagation along n and (n, u, v) cyclic:
 * @verbatim
      E_u(p)     += dt / (eps dn) * H_v,inc(p - 1/2)      H_v(p - 1/2) += dt / (mu0 dn) * E_u,inc(p)
      E_v(p)     -= dt / (eps dn) * H_u,inc(p - 1/2)      H_u(p - 1/2) -= dt / (mu0 dn) * E_v,inc(p)
   @endverbatim
 * The incident fields are the mode profile times the signal, with the delay of half a cell at
 * n_eff between the E and H planes. The sign of H is chosen so that the mode's power flows along +n.
 * See docs/sources/ModeSource.md.
 */
#ifndef FDTD_SOURCES_MODE_SOURCE_H
#define FDTD_SOURCES_MODE_SOURCE_H

#include "sources/mode_spec.h"
#include "sources/signal.h"

namespace fdtd::sources
{
extern ModeSpec mode_input; ///< [mode], set by apply_sim_config

/// Reads the MPB files and builds this rank's profiles; fills the full material maps, PML included
void setup_mode_source();
void mode_source_h(const Signal &signal, double t); ///< after the H update; t: the step's time
void mode_source_e(const Signal &signal, double t); ///< after the E update
} // namespace fdtd::sources

#endif
