/**
 * @file sources/state.h
 * @brief State of the sources module
 *
 * Moved out of the former shared namespace fdtd::state: read and written by this module, allocated or set by the
 * run set-up, app/setup.cc, and the configuration, applyConfig.cc. Defined in src/sources/state.cc.
 */
#ifndef FDTD_SOURCES_STATE_H
#define FDTD_SOURCES_STATE_H

#include "core/carray.h"
#include "sources/types.h"
#include "core/constants.h"
#include "materials/material_field.h"
#include "materials/codes.h"
#include <complex>
#include <mpi.h>
#include <string>
#include <vector>
#include "core/limits.h"

namespace fdtd::sources
{
extern double chx;                                            ///< TF/SF H correction coefficient DT / (MU0 * DX) used at the TF/SF boundary.
extern double chy;                                            ///< TF/SF H correction coefficient DT / (MU0 * DY) used at the TF/SF boundary.
extern double chz;                                            ///< TF/SF H correction coefficient DT / (MU0 * DZ) used at the TF/SF boundary.
extern int jpw_start_uni;                                     ///< 1D grid index sampled to drive the E side of the TF/SF boundary.
extern int jpw_source_uni;                                    ///< Index in the 1D auxiliary grid where the plane-wave source is injected.
extern int jpw_end_uni;                                     ///< 1D grid index sampled at the far TF/SF boundary.
extern int nr_cell_par_uni;                                 ///< 1D grid: first cell of the far-side layer, second medium.
extern int j_min_cpml_uni;                                  ///< 1D grid: first index of the middle region between the two CPMLs.
extern int j_max_cpml_uni;                                  ///< 1D grid: last index of the middle region between the two CPMLs.
extern int nr_cell_uni;                                     ///< 1D grid: last cell index (nr_cells_uni - 1).
extern int ipw_near_h;                                           ///< Local x index of the near TF/SF plane, H side.
extern int jpw_near_h;                                           ///< Local y index of the near TF/SF plane, H side.
extern int kpw_near_h;                                           ///< Local z index of the near TF/SF plane, H side.
extern int ipw_min_mod;                                     ///< First x index of the plane-wave injection line, adjusted for the field component.
extern int kpw_min_mod;                                     ///< First z index of the plane-wave injection line, adjusted for the field component.
extern double polarization_angle;                             ///< Plane-wave polarisation angle in the xz plane in rad.
extern std::vector<PlaneWave> plane_wave_input;                  ///< Plane waves as given in the configuration.
extern std::vector<DipoleInput> dipole_input;               ///< Point dipoles as given in the configuration.
extern bool plane_wave_dispersion;                          ///< Correct the plane-wave source for numerical dispersion.
extern int n_start;                                         ///< First time step included in the source-signal DFT.
extern std::vector<double> spectrum_omega;                    ///< Angular frequencies of the spectrum points in rad/s.
extern double fundamental_step_frequency;                     ///< Frequency step of the fundamental band in Hz.
extern bool xpol;                                           ///< Source polarised along x.
extern bool zpol;                                           ///< Source polarised along z.
extern std::string pol;                                          ///< Polarisation label used in output file names.
extern double input_power;                                    ///< Incident fluence integrated over the spectrum in J/m^2.
extern double analysis_omega[fdtd::max_frequencies];                      ///< Analysed angular frequencies in rad/s.
extern double amplitude_normalization[fdtd::max_frequencies];    ///< |DFT of the source signal| per frequency, used to normalise spectra.
extern std::complex<double> dft_normalization[fdtd::max_frequencies];         ///< Complex DFT of the source signal per frequency.
} // namespace fdtd::sources

#endif
