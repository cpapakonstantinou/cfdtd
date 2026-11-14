/**
 * @file sources/superposition.h
 * @brief Plane-wave superposition source, wave_superposition, Gaussian beam, on the TF/SF faces
 *
 * Each wave of the superposition, ::PlaneWave, inc/Sources.h: direction, E and H amplitudes, pulse or
 * its derivative, grid-dispersion factor, is injected analytically on the six TF/SF faces: every
 * tangential component on the face gets amplitude * signal(t - delay), the delay being the travel
 * time from the focus to that component's Yee position in the first medium,
 * @verbatim
      delay = (d . (r - r_focus)) n_first / (c0 dispersion),   summed over x, y, z in that order
   @endverbatim
 * H is corrected with the incident E on the face, E with the incident H half a cell outside it.
 */
#ifndef FDTD_SOURCES_SUPERPOSITION_H
#define FDTD_SOURCES_SUPERPOSITION_H

#include "sources/types.h"
#include "sources/signal.h"
#include <vector>

namespace fdtd::sources
{
/// H corrections of every wave at time t on this rank's TF/SF faces
void superposition_h(const Signal &signal, const WaveFrame &frame, double t, const std::vector<PlaneWave> &waves);
/// E corrections of every wave at time t on this rank's TF/SF faces
void superposition_e(const Signal &signal, const WaveFrame &frame, double t, const std::vector<PlaneWave> &waves);
} // namespace fdtd::sources

#endif
