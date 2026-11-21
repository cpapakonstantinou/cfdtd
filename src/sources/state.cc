/**
 * @file sources/state.cc
 * @brief Definitions of the sources module state, sources/state.h
 */
#include "sources/state.h"
#include "sources/spec.h"
#include "core/limits.h"

namespace fdtd::sources
{
double chx;
double chy;
double chz;
int jpw_start_uni;
int jpw_source_uni;
int jpw_end_uni;
int nr_cell_par_uni;
int j_min_cpml_uni;
int j_max_cpml_uni;
int nr_cell_uni;
int ipw_near_h;
int jpw_near_h;
int kpw_near_h;
int ipw_min_mod;
int kpw_min_mod;
double polarization_angle;
std::vector<PlaneWave> plane_wave_input;
std::vector<DipoleInput> dipole_input;
bool plane_wave_dispersion = true;
int n_start;
std::vector<double> spectrum_omega;
double fundamental_step_frequency;
bool xpol;
bool zpol;
double input_power;
double analysis_omega[fdtd::max_frequencies];
double amplitude_normalization[fdtd::max_frequencies];
std::complex<double> dft_normalization[fdtd::max_frequencies];
} // namespace fdtd::sources

fdtd::sources::Spec fdtd::sources::the_spec;
