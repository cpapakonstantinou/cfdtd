/** \file   core/quadrature.h
 *  \brief  Gauss-Legendre quadrature: nodes and weights on [-1, 1], and a panel rule on [a, b].
 */
#ifndef FDTD_CORE_QUADRATURE_H
#define FDTD_CORE_QUADRATURE_H

#include <vector>

namespace fdtd::core
{
/// Gauss-Legendre nodes, ascending, and weights on [-1, 1], n >= 1; exact for polynomials of degree 2n - 1
void gauss_legendre(int n, std::vector<double> &x, std::vector<double> &w);

/// nodes and weights of `panels` equal Gauss-Legendre panels of `n` nodes each on [a, b]
void gauss_panels(double a, double b, int panels, int n, std::vector<double> &x, std::vector<double> &w);
} // namespace fdtd::core

#endif // FDTD_CORE_QUADRATURE_H
