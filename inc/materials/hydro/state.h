/**
 * @file materials/hydro/state.h
 * @brief State of the electron-gas models: hydrodynamic densities, currents, margin grid, coefficients
 *
 * Moved out of the former shared namespace fdtd::state: read and written by this module, allocated or set by the
 * run set-up, app/setup.cc, and the configuration, applyConfig.cc. Defined in src/materials/hydro/state.cc.
 */
#ifndef FDTD_MATERIALS_HYDRO_STATE_H
#define FDTD_MATERIALS_HYDRO_STATE_H

#include "core/carray.h"
#include "core/constants.h"
#include "materials/material_field.h"
#include "materials/codes.h"
#include <complex>
#include <mpi.h>
#include <string>
#include <vector>
#include "materials/dispersive_state.h"

namespace fdtd::materials::hydro
{
inline ctensor<double> &ex_2_nl = fdtd::materials::the_dispersive_state.c[0].e_padded_nl; ///< Copy of the harmonic Ex on the margin grid.
inline ctensor<double> &ey_2_nl = fdtd::materials::the_dispersive_state.c[1].e_padded_nl; ///< Copy of the harmonic Ey on the margin grid.
inline ctensor<double> &ez_2_nl = fdtd::materials::the_dispersive_state.c[2].e_padded_nl; ///< Copy of the harmonic Ez on the margin grid.
inline material_field<double, materials::nr_disp_mat_d> &jx_prev = fdtd::materials::the_dispersive_state.c[0].j_prev;        ///< Jx at the previous step.
inline material_field<double, materials::nr_disp_mat_d> &jy_prev = fdtd::materials::the_dispersive_state.c[1].j_prev;        ///< Jy at the previous step.
inline material_field<double, materials::nr_disp_mat_d> &jz_prev = fdtd::materials::the_dispersive_state.c[2].j_prev;        ///< Jz at the previous step.
inline material_field<double, materials::nr_disp_mat_d> &jx_extra = fdtd::materials::the_dispersive_state.c[0].j_extra;       ///< Hydrodynamic source term of Jx, convection, Lorentz force, pressure.
inline material_field<double, materials::nr_disp_mat_d> &jy_extra = fdtd::materials::the_dispersive_state.c[1].j_extra;       ///< Hydrodynamic source term of Jy, convection, Lorentz force, pressure.
inline material_field<double, materials::nr_disp_mat_d> &jz_extra = fdtd::materials::the_dispersive_state.c[2].j_extra;       ///< Hydrodynamic source term of Jz, convection, Lorentz force, pressure.
inline ctensor<double> &nd_x = fdtd::materials::the_dispersive_state.c[0].electron_density; ///< Induced electron density from the total fields, margin grid, m^-3.
inline ctensor<double> &nd_y = fdtd::materials::the_dispersive_state.c[1].electron_density; ///< Induced electron density from the harmonic fields, margin grid, m^-3.
inline ctensor<double> &nd_z = fdtd::materials::the_dispersive_state.c[2].electron_density; ///< Induced electron density, second estimate, margin grid, m^-3.
extern ctensor<double> nd_z_prev;                            ///< Induced electron density at the previous step.
extern bool with_pressure;                                  ///< Include the electron pressure, Thomas-Fermi, term.
extern bool with_mag_field;                                 ///< Include the magnetic Lorentz force on the electrons.
extern bool with_convection;                                ///< Include the convective (v.grad)v term.
inline ctensor<double> &ex_2 = fdtd::materials::the_dispersive_state.c[0].e_padded; ///< Copy of Ex on the margin grid, cells+2, used for the density stencil.
inline ctensor<double> &ey_2 = fdtd::materials::the_dispersive_state.c[1].e_padded; ///< Copy of Ey on the margin grid, cells+2, used for the density stencil.
inline ctensor<double> &ez_2 = fdtd::materials::the_dispersive_state.c[2].e_padded; ///< Copy of Ez on the margin grid, cells+2, used for the density stencil.
inline ctensor<double> &epsinf_x = fdtd::materials::the_dispersive_state.c[0].eps_inf; ///< eps_inf at Ex positions on the margin grid.
inline ctensor<double> &epsinf_y = fdtd::materials::the_dispersive_state.c[1].eps_inf; ///< eps_inf at Ey positions on the margin grid.
inline ctensor<double> &epsinf_z = fdtd::materials::the_dispersive_state.c[2].eps_inf; ///< eps_inf at Ez positions on the margin grid.
inline ctensor<double> &grad_div_x = fdtd::materials::the_dispersive_state.c[0].grad_div; ///< Nonlocal term grad(div P), x component.
inline ctensor<double> &grad_div_y = fdtd::materials::the_dispersive_state.c[1].grad_div; ///< Nonlocal term grad(div P), y component.
inline ctensor<double> &grad_div_z = fdtd::materials::the_dispersive_state.c[2].grad_div; ///< Nonlocal term grad(div P), z component.
inline ctensor<double> &grad_div_x_n = fdtd::materials::the_dispersive_state.c[0].grad_div_n; ///< grad(div P) x at the current step, GNOR history.
inline ctensor<double> &grad_div_y_n = fdtd::materials::the_dispersive_state.c[1].grad_div_n; ///< grad(div P) y at the current step, GNOR history.
inline ctensor<double> &grad_div_z_n = fdtd::materials::the_dispersive_state.c[2].grad_div_n; ///< grad(div P) z at the current step, GNOR history.
inline ctensor<double> &grad_div_x_n_1 = fdtd::materials::the_dispersive_state.c[0].grad_div_n_1; ///< grad(div P) x one step back, GNOR history.
inline ctensor<double> &grad_div_y_n_1 = fdtd::materials::the_dispersive_state.c[1].grad_div_n_1; ///< grad(div P) y one step back, GNOR history.
inline ctensor<double> &grad_div_z_n_1 = fdtd::materials::the_dispersive_state.c[2].grad_div_n_1; ///< grad(div P) z one step back, GNOR history.
inline ctensor<double> &px_nl = fdtd::materials::the_dispersive_state.c[0].p_nonlocal; ///< Nonlocal, hydrodynamic, polarisation Px.
inline ctensor<double> &py_nl = fdtd::materials::the_dispersive_state.c[1].p_nonlocal; ///< Nonlocal, hydrodynamic, polarisation Py.
inline ctensor<double> &pz_nl = fdtd::materials::the_dispersive_state.c[2].p_nonlocal; ///< Nonlocal, hydrodynamic, polarisation Pz.
inline ctensor<double> &px_nl_2 = fdtd::materials::the_dispersive_state.c[0].p_nonlocal_2; ///< Second nonlocal polarisation Px, Laplacian form.
inline ctensor<double> &py_nl_2 = fdtd::materials::the_dispersive_state.c[1].p_nonlocal_2; ///< Second nonlocal polarisation Py, Laplacian form.
inline ctensor<double> &pz_nl_2 = fdtd::materials::the_dispersive_state.c[2].p_nonlocal_2; ///< Second nonlocal polarisation Pz, Laplacian form.
inline ctensor<double> &px_cp_2 = fdtd::materials::the_dispersive_state.c[0].pc_sum; ///< Sum of the critical-point Px on the margin grid.
inline ctensor<double> &py_cp_2 = fdtd::materials::the_dispersive_state.c[1].pc_sum; ///< Sum of the critical-point Py on the margin grid.
inline ctensor<double> &pz_cp_2 = fdtd::materials::the_dispersive_state.c[2].pc_sum; ///< Sum of the critical-point Pz on the margin grid.
inline ctensor<double> &px_cp_2_nl = fdtd::materials::the_dispersive_state.c[0].pc_sum_nl; ///< Sum of the harmonic critical-point Px on the margin grid.
inline ctensor<double> &py_cp_2_nl = fdtd::materials::the_dispersive_state.c[1].pc_sum_nl; ///< Sum of the harmonic critical-point Py on the margin grid.
inline ctensor<double> &pz_cp_2_nl = fdtd::materials::the_dispersive_state.c[2].pc_sum_nl; ///< Sum of the harmonic critical-point Pz on the margin grid.
inline double (&non_loc)[materials::nr_disp_mat_d] = fdtd::materials::the_model_parameters.non_loc; ///< Nonlocal, hydrodynamic, coupling beta^2 term per material.
inline double (&beta_nl)[materials::nr_disp_mat_d] = fdtd::materials::the_model_parameters.beta_nl; ///< Hydrodynamic beta in m/s per material.
inline double (&diff_coeff)[materials::nr_disp_mat_d] = fdtd::materials::the_model_parameters.diff_coeff; ///< GNOR diffusion coefficient per material in m^2/s.
inline double (&gnor_1)[materials::nr_disp_mat_d] = fdtd::materials::the_model_parameters.gnor_1; ///< GNOR coefficient on grad(div P) at step n.
inline double (&gnor_2)[materials::nr_disp_mat_d] = fdtd::materials::the_model_parameters.gnor_2; ///< GNOR coefficient on the grad(div P) time difference.
inline double (&equilibrium_electron_density)[materials::nr_disp_mat_d] = fdtd::materials::the_model_parameters.equilibrium_electron_density; ///< Equilibrium free-electron density per material in m^-3.
inline double (&effective_electron_mass)[materials::nr_disp_mat_d] = fdtd::materials::the_model_parameters.effective_electron_mass; ///< Effective electron mass per material (kg).
inline double (&powerterm)[materials::NRMAXMATERIALS] = fdtd::materials::the_model_parameters.powerterm; ///< n0^(2/3) per material, for the pressure term.
inline double (&inverse_edens)[materials::NRMAXMATERIALS] = fdtd::materials::the_model_parameters.inverse_edens; ///< 1 / n0 per material.
inline double (&inverse_chargedens)[materials::NRMAXMATERIALS] = fdtd::materials::the_model_parameters.inverse_chargedens; ///< 1 / (e n0) per material.
extern bool use_backwards_differences;                      ///< Use backward differences for the density gradient.
extern int nonlocal_barrier_int;                            ///< Global y index below which the nonlocal term is applied.
} // namespace fdtd::materials::hydro

#endif
