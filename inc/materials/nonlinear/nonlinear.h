/**
 * @file materials/nonlinear/nonlinear.h
 * @brief Second- and third-order nonlinearities
 *
 * @verbatim
      update_polarisation()   P_NL approach: P += P_NL(E) on the nonlinear cells of the dispersive
                              nonlinear materials, before the E scan, tensor.cc
      update_field()          E approach: E += the change of P_NL over the step / eps, after the E
                              scan, tensor.cc
      update_exact_shg/thg    dispersionless isotropic chi2 / chi3 dielectrics: the E update of one
                              cell solved in closed form, exact.cc
      polarisation_term       the chi contraction of one component, below
   @endverbatim
 */
#ifndef FDTD_MATERIALS_NONLINEAR_NONLINEAR_H
#define FDTD_MATERIALS_NONLINEAR_NONLINEAR_H

namespace fdtd::materials::nonlinear
{
void update_polarisation();
void update_field();
/// E_a at (i, j, k) of the dielectric material m from D = eps E + chi E^2, SHG, or chi E^3, THG;
/// rotdiff: the curl term of the step
void update_exact_shg(int a, double rotdiff, int m, int i, int j, int k);
void update_exact_thg(int a, double rotdiff, int m, int i, int j, int k);

/**
 * @brief Nonlinear polarisation of component a (0 = x, 1 = y, 2 = z) at (i, j, k), or its change over one step
 * @verbatim
      chi2, SHG, Voigt row 6a:    P_a = eps0 sum_q d[6a + q] (EE)_q,   (EE) = (xx, yy, zz, 2yz, 2xz, 2xy)
      chi3, THG, row 7a:          P_a = eps0 E_a (c[7a] E_a^2 + (c[7a+1..3]) E_b^2 + (c[7a+4..6]) E_c^2)
                                   b < c: the other two axes
      change = true:               every product E_p E_q becomes E_p E_q - E_p,prev E_q,prev, field approach
   @endverbatim
 * Staircased cells whose three components share one material use the full tensor. Otherwise only the
 * terms of components whose material matches component a's are kept: a component of another material
 * does not belong to this medium. With SHG and THG both set, the SHG result is returned, as before.
 *
 * @param mm       dispersive model index of the nonlinear material, chi tables chi_2[mm], chi_3[mm]
 * @param mx, my, mz  materials of the three components at the cell (matind_x/y/z)
 */
double polarisation_term(int mm, int i, int j, int k, int a, int mx, int my, int mz, bool change);
} // namespace fdtd::materials::nonlinear

#endif
