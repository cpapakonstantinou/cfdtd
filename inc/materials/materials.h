/**
 * @file materials/materials.h
 * @brief The material table: parameters of every material code, its optical constants at a
 *        frequency, and the ADE update coefficients derived from them
 *
 * @verbatim
      set_parameters()                       fills the per-code model parameters, Drude, Lorentz,
                                             critical points, plasma frequencies, ..., from the config
      refractive_index(code, w)              n(w) = sqrt(eps_r(w))
      complex_permittivity(code, w)          eps_r(w) of the code's model
      compute_ade_coefficients()             the per-code coefficients of the E and P updates for DT
   @endverbatim
 * Material codes and their ranges: materials/codes.h.
 */
#ifndef FDTD_MATERIALS_MATERIALS_H
#define FDTD_MATERIALS_MATERIALS_H

#include <string>

#include <complex>

namespace fdtd::materials
{
struct Inputs;
void set_parameters(const Inputs &in);
std::complex<double> refractive_index(int code, double omega);
std::complex<double> complex_permittivity(int code, double omega);
void compute_ade_coefficients();
/// the static magnetisation field of the magnetic Lorentz model, from its input file, magnetisation.cc
void read_magnetisation_field(const std::string &input_dir); ///< Magnetization_{x,y,z}.txt in input_dir
} // namespace fdtd::materials

#endif
