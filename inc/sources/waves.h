/**
 * @file sources/waves.h
 * @brief The plane-wave list of the superposition and Gaussian-beam sources
 */
#ifndef FDTD_SOURCES_WAVES_H
#define FDTD_SOURCES_WAVES_H

#include "sources/types.h"
#include <vector>

namespace fdtd::sources
{
/**
 * @brief Builds the waves of the superposition source: the [wave.*] sections, or the angular
 *        decomposition of the Gaussian beam, or one oblique wave when it has no angular steps
 *
 * Also sets what the face corrections read: the focus, the TF/SF box positions (x/y/z_tfsf_n/f, m)
 * and the beam's time offset. Nothing to do, true, empty list, unless the source is a superposition
 * or a Gaussian beam.
 * @return false when the configured waves are invalid
 */
bool build_waves(std::vector<PlaneWave> &waves, WaveFrame &frame);
} // namespace fdtd::sources

#endif
