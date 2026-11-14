/**
 * @file monitors/state.h
 * @brief State of the monitors module
 *
 * Moved out of the former shared namespace fdtd::state: read and written by this module, allocated or set by the
 * run set-up, app/setup.cc, and the configuration, applyConfig.cc. Defined in src/monitors/state.cc.
 */
#ifndef FDTD_MONITORS_STATE_H
#define FDTD_MONITORS_STATE_H

#include "core/carray.h"
#include "core/constants.h"
#include "materials/material_field.h"
#include "materials/codes.h"
#include <complex>
#include <mpi.h>
#include <string>
#include <vector>
#include "core/limits.h"

namespace fdtd::monitors
{
extern double max_source_frequency;                                   ///< Highest source frequency in Hz.
extern bool dft_plane_xz;                                         ///< Record DFTs on the xz output plane.
extern bool dft_xz_1;                                       ///< Record DFTs on xz output plane 1.
extern bool dft_xz_2;                                       ///< Record DFTs on xz output plane 2.
extern bool dft_plane_yz;                                         ///< Record DFTs on the yz output plane.
extern bool dft_plane_xy;                                         ///< Record DFTs on the xy output plane.
extern bool dft_h;                                          ///< Also record H on the output planes.
extern bool phase_normalization;                            ///< Remove the source phase when normalising DFTs.
extern double harmonic_step_frequency;                        ///< Frequency step of the harmonic band in Hz.
extern FILE *f_art;                                         ///< Output file of absorption/reflection/transmission.
extern double transm_nl_power;                                ///< Harmonic energy scattered or transmitted in J.
extern double refl_nl_power;                                  ///< Harmonic energy reflected or absorbed in J.
extern double abs_sigma_loc[fdtd::max_frequencies];             ///< This rank's part of abs_sigma_spectrum.
extern double abs_sigma_spectrum[fdtd::max_frequencies];                  ///< Absorbed power per frequency from sigma |E|^2 over the volume.
extern double abs_sigma_2[fdtd::max_frequencies];                ///< Second absorption estimate per frequency from sigma |E|^2.
extern double area_loc;                                       ///< This rank's part of Area.
extern double q_normalisation_area;                                           ///< Area that normalises the Q coefficients in m^2.
extern bool q_coeff;                                        ///< Compute scattering/absorption cross sections, Q coefficients.
extern bool single_photon;                                  ///< Normalise the harmonic output per incident photon.
extern std::string str_2;                                        ///< Scratch: run label used in ART file names.
extern long int startindices_nf[3];                         ///< Global index of this rank's first near-field cell (x y z).
extern long int globalsizes_nf[3];                          ///< Global near-field block size, cells, x y z.
extern long int localsize_nf[3];                            ///< Size of this rank's near-field block, cells, x y z.
extern MPI_Comm output_comm_sca_nea_x;                      ///< Communicator of the ranks on the near x face of the scattering box.
extern MPI_Comm output_comm_sca_nea_y;                      ///< Communicator of the ranks on the near y face of the scattering box.
extern MPI_Comm output_comm_sca_nea_z;                      ///< Communicator of the ranks on the near z face of the scattering box.
extern MPI_Comm output_comm_sca_far_x;                      ///< Communicator of the ranks on the far x face of the scattering box.
extern MPI_Comm output_comm_sca_far_y;                      ///< Communicator of the ranks on the far y face of the scattering box.
extern MPI_Comm output_comm_sca_far_z;                      ///< Communicator of the ranks on the far z face of the scattering box.
extern int sphere_rad_multipole;                            ///< Radius of the multipole integration sphere in cells.
extern int ntheta_multipole;                                ///< Number of polar angles on the multipole sphere.
extern int nphi_multipole;                                  ///< Number of azimuthal angles on the multipole sphere.
extern bool exact_multipoles;                               ///< Compute multipoles from the exact, volume current, formulas.
extern bool lwa_multipoles;                                 ///< Compute the long-wavelength (kr -> 0) multipoles from the same volume current.
inline bool multipoles_wanted() { return exact_multipoles || lwa_multipoles; } ///< either multipole set: needs the volume E-field DFT
extern MPI_Comm output_comm_plot_xz_1;                      ///< Communicator of the ranks that cut xz output plane 1.
extern int count_art;                                       ///< Counter of ART evaluations; output is written every 10th.
extern double geometric_cross_section;                        ///< Geometric cross section of the drawn shape in m^2; -1 if the shape does not set it.
extern ctensor<double> energy_density;                        ///< Electromagnetic energy density per cell in J/m^3.
} // namespace fdtd::monitors

#endif
