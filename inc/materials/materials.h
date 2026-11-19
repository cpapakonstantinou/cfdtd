/**
 * @file materials/materials.h
 * @brief The materials of a run: optical constants of a slot at a frequency, and the ADE update
 *        coefficients derived from the parameters
 *
 * @verbatim
      load_slots(definitions)                the slots and their parameters, materials/slots.h
      refractive_index(slot, w)              n(w) = sqrt(eps_r(w))
      complex_permittivity(slot, w)          eps_r(w) of the slot's model
      compute_ade_coefficients()             the coefficients of the E and P updates for the time step
   @endverbatim
 */
#ifndef FDTD_MATERIALS_MATERIALS_H
#define FDTD_MATERIALS_MATERIALS_H

#include <string>

#include <complex>

namespace fdtd::materials
{
std::complex<double> refractive_index(int s, double omega);
std::complex<double> complex_permittivity(int s, double omega);
void compute_ade_coefficients();
/// the static magnetisation field of the magnetic Lorentz model, from its input file, magnetisation.cc
void read_magnetisation_field(const std::string &input_dir); ///< Magnetization_{x,y,z}.txt in input_dir
} // namespace fdtd::materials

#endif
