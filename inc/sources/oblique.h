/**
 * @file sources/oblique.h
 * @brief Plane waves through a layered background on the TF/SF faces
 *
 * When a superposition source meets a layered background, middle or last medium differs from the
 * first, the incident field is taken from 1D simulations of each wave through the layers
 * OneDimensionalSimulator: E and H at any point and time. The face corrections are those of
 * sources/superposition.h, Taflove ch. 5, with the E coefficient of the layer at the cell.
 */
#ifndef FDTD_SOURCES_OBLIQUE_H
#define FDTD_SOURCES_OBLIQUE_H

class OneDimensionalSimulator;

namespace fdtd::sources
{
/// H corrections at time t on this rank's TF/SF faces
void oblique_h(OneDimensionalSimulator &ods, double t);
/// E corrections at time t on this rank's TF/SF faces
void oblique_e(OneDimensionalSimulator &ods, double t);
} // namespace fdtd::sources

#endif
