/** \file   core/constants.h
 *  \brief  Physical constants, SI units, and solver-wide fixed values.
 */
#ifndef FDTD_CONSTANTS_H
#define FDTD_CONSTANTS_H

#include <cmath>
#include <complex>

namespace fdtd::constants
{
inline constexpr double EP0 = 8.85418781762E-12; ///< vacuum permittivity in F/m
inline constexpr double MU0 = 1.2566370614E-06;  ///< vacuum permeability (H/m)
inline const double C0 = 1 / (std::sqrt(EP0 * MU0)); ///< speed of light in m/s, from EP0 and MU0
inline const double INVVACIMP = std::sqrt(EP0 / MU0); ///< 1 / vacuum impedance (S)
inline constexpr double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
inline constexpr double PI2 = 2. * PI;
inline constexpr double PI3 = 3. * PI;
inline constexpr double E0 = 1.60217662E-19;   ///< elementary charge in C
inline constexpr double ME = 9.10938356E-31;   ///< electron mass (kg)
inline constexpr double HBAR = 1.0545718E-34;  ///< reduced Planck constant in J s
inline constexpr std::complex<double> I_UNIT{0., 1.}; ///< imaginary unit
inline constexpr double DEFAULT_TEMPERATURE = 300.0; ///< initial electron and lattice temperature in K
inline constexpr int ROOTPROCESS = 0;                ///< rank that reads input and writes shared output
} // namespace fdtd::constants

#endif // FDTD_CONSTANTS_H
