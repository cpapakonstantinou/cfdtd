/**
 * @file materials/nonlinear/state.cc
 * @brief Definitions of the nonlinear module state, nonlinear/state.h
 */
#include "materials/nonlinear/state.h"

namespace fdtd::materials::nonlinear
{
using namespace fdtd::materials;

double chi_2[nr_disp_mat_d][18];
double chi_2_a[18];
double chi_2_f[18];
double chi_2_m[18];
double chi_2_s[18];
double chi_3[nr_disp_mat_d][21];
double chi_3_a[21];
double chi_3_f[21];
double chi_3_m[21];
double chi_3_s[21];
double chi2_dl[DispMateOffsetD][18];
double chi3_dl[DispMateOffsetD][21];
bool iso_chi_2_dl[DispMateOffsetD];
bool iso_chi_3_dl[DispMateOffsetD];
bool nl[nr_disp_mat_d];
bool nl_dl[DispMateOffsetD];
bool is_chi_2_a_d;
bool is_chi_2_a_o;
bool is_chi_3_a_d;
bool is_chi_3_a_o;
} // namespace fdtd::materials::nonlinear
