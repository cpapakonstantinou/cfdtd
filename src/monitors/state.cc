/**
 * @file monitors/state.cc
 * @brief Definitions of the monitors module state, monitors/state.h
 */
#include "monitors/state.h"
#include "monitors/spec.h"
#include "core/limits.h"

namespace fdtd::monitors
{
double max_source_frequency;
bool dft_plane_xz;
bool dft_xz_1;
bool dft_xz_2;
bool dft_plane_yz;
bool dft_plane_xy;
bool dft_h;
bool phase_normalization;
double harmonic_step_frequency;
FILE *f_art = NULL;
double transm_nl_power;
double refl_nl_power;
double abs_sigma_loc[fdtd::max_frequencies];
double abs_sigma_spectrum[fdtd::max_frequencies];
double abs_sigma_2[fdtd::max_frequencies];
double area_loc;
double q_normalisation_area;
bool q_coeff;
bool single_photon;
std::string str_2;
long int startindices_nf[3];
long int globalsizes_nf[3];
long int localsize_nf[3];
MPI_Comm output_comm_sca_nea_x;
MPI_Comm output_comm_sca_nea_y;
MPI_Comm output_comm_sca_nea_z;
MPI_Comm output_comm_sca_far_x;
MPI_Comm output_comm_sca_far_y;
MPI_Comm output_comm_sca_far_z;
int sphere_rad_multipole;
int ntheta_multipole;
int nphi_multipole;
bool exact_multipoles = false;
bool lwa_multipoles = false;
MPI_Comm output_comm_plot_xz_1;
int count_art;
double geometric_cross_section = -1;
ctensor<double> energy_density;
} // namespace fdtd::monitors

fdtd::monitors::Spec fdtd::monitors::the_spec;
