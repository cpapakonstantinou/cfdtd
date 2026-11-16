/**
 * @file materials/coefficients.h
 * @brief The coefficient table of each material model
 *
 * A model type M declares `struct Coeffs`, its precomputed update coefficients for one material, and
 * `static Coeffs compute(int mm, int code, double dt)`, the expressions, from the material parameters
 * and the time step. `compute_ade_coefficients()` fills `coeff_table<M>[mm]` for every dispersive
 * material index mm of M's code range; the update and every other reader take `coeffs<M>(mm)`.
 * Entries of indices outside M's range stay zero.
 */
#ifndef FDTD_MATERIALS_COEFFICIENTS_H
#define FDTD_MATERIALS_COEFFICIENTS_H

#include "materials/codes.h"
#include <array>

namespace fdtd::materials
{
/**
 * @brief Critical-point, or second-order, poles and the E update of the models built on them
 * @verbatim
      Pn_new = b0_n Pn_prev + b1_n Pn + b2_n E_prev + b3_n E + b4_n E_new
      E_new  = c0 curl + c1 (...) + c2 E_prev + c3 E, the bracket differs per model
   @endverbatim
 */
struct Poles
{
	double b0[NrMaxCP];
	double b1[NrMaxCP];
	double b2[NrMaxCP];
	double b3[NrMaxCP];
	double b4[NrMaxCP];
	double c0;
	double c1;
	double c2;
	double c3;
};

/// free-electron current J_new = nl1 J + nl2 E (+ nl3 n1 E, nl4 and nl5 terms of the hydrodynamic model)
struct ElectronGas
{
	double nl1;
	double nl2;
	double nl3;
	double nl4;
	double nl5;
};

/// coefficients of model M per dispersive material index, code - DispMateOffsetD
template <class M>
inline std::array<typename M::Coeffs, nr_disp_mat_d> coeff_table{};

/// the coefficients of material index mm under model M
template <class M>
inline const typename M::Coeffs &coeffs(int mm) { return coeff_table<M>[mm]; }
} // namespace fdtd::materials

#endif
