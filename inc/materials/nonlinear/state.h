/**
 * @file materials/nonlinear/state.h
 * @brief State of the nonlinear module
 *
 * Moved out of the former shared namespace fdtd::state: read and written by this module, allocated or set by the
 * run set-up, app/setup.cc, and the configuration, applyConfig.cc. Defined in src/materials/nonlinear/state.cc.
 */
#ifndef FDTD_MATERIALS_NONLINEAR_STATE_H
#define FDTD_MATERIALS_NONLINEAR_STATE_H

#include "core/carray.h"
#include "core/constants.h"
#include "materials/material_field.h"
#include "materials/codes.h"
#include <complex>
#include <mpi.h>
#include <string>
#include <vector>

namespace fdtd::materials::nonlinear
{
extern double chi_2[materials::nr_disp_mat_d][18];                        ///< chi(2) tensor, 18 Voigt entries, per dispersive material.
extern double chi_2_a[18];                                    ///< chi(2) tensor of the antenna as read from the configuration.
extern double chi_2_f[18];                                    ///< chi(2) tensor of the first medium as read from the configuration.
extern double chi_2_m[18];                                    ///< chi(2) tensor of the middle layer as read from the configuration.
extern double chi_2_s[18];                                    ///< chi(2) tensor of the second medium as read from the configuration.
extern double chi_3[materials::nr_disp_mat_d][21];                       ///< chi(3) tensor, 21 entries, per dispersive material.
extern double chi_3_a[21];                                    ///< chi(3) tensor of the antenna as read from the configuration.
extern double chi_3_f[21];                                    ///< chi(3) tensor of the first medium as read from the configuration.
extern double chi_3_m[21];                                    ///< chi(3) tensor of the middle layer as read from the configuration.
extern double chi_3_s[21];                                    ///< chi(3) tensor of the second medium as read from the configuration.
extern double chi2_dl[materials::DispMateOffsetD][18];                   ///< chi(2) tensor, 18 Voigt entries, per dielectric material code.
extern double chi3_dl[materials::DispMateOffsetD][21];                   ///< chi(3) tensor, 21 entries, per dielectric material code.
extern bool iso_chi_2_dl[materials::DispMateOffsetD];                  ///< chi(2) of each dielectric material is isotropic, scalar.
extern bool iso_chi_3_dl[materials::DispMateOffsetD];                  ///< chi(3) of each dielectric material is isotropic, scalar.
extern bool nl[materials::nr_disp_mat_d];                              ///< Dispersive material is nonlinear.
extern bool nl_dl[materials::DispMateOffsetD];                         ///< Dielectric material is nonlinear.
extern bool is_chi_2_a_d;                                   ///< chi(2) of the antenna has a nonzero diagonal entry.
extern bool is_chi_2_a_o;                                   ///< chi(2) of the antenna has a nonzero off-diagonal entry.
extern bool is_chi_3_a_d;                                   ///< chi(3) of the antenna has a nonzero diagonal entry.
extern bool is_chi_3_a_o;                                   ///< chi(3) of the antenna has a nonzero off-diagonal entry.
} // namespace fdtd::materials::nonlinear

#endif
