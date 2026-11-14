/**
 * @file materials/nonlinear/exact.cc
 * @brief Exact E update of dispersionless, isotropic chi2 / chi3 dielectrics, SHG and THG
 *
 * For a dispersionless isotropic medium the update D = eps0 (eps E + chi E^n) is a quadratic, SHG, or
 * cubic, THG, equation in the new E, solved in closed form. One solver per order serves Ex, Ey and Ez:
 * the component's field array and its diagonal chi element, Voigt xx, yy, zz, are the arguments.
 * The chi tables are the dielectric ones, chi2_dl, chi3_dl, index = material code.
 */
#include "materials/nonlinear/nonlinear.h"
#include "fields/fields.h"
#include "core/views.h"
#include "materials/nonlinear/state.h"
#include "materials/state.h"
#include "materials/spec.h"
using fdtd::constants::EP0;

namespace
{
} // namespace

namespace
{
/**
 * @brief SHG: solve chi E^2 + eps E - D = 0 for the new E
 * @verbatim
      D = eps E + chi E^2 + dt curl / eps0
      |4 chi D| < 1e-8 eps^2:    E = D / eps, linear limit
      eps^2 + 4 chi D >= 0:      E = (-eps + sqrt(eps^2 + 4 chi D)) / (2 chi)
      otherwise:                 E = D / (eps + chi E_old), no real root
   @endverbatim
 */
void exact_shg(view3<double> e, int m, const double *chi2, double rotdiff, int i, int j, int k)
{
	const auto &mp = fdtd::materials::model_parameters();
	const double chi = chi2[0];
	const double d = mp.epre[m] * e[i][j][k] + chi * e[i][j][k] * e[i][j][k] + fdtd::materials::model_spec().dt * rotdiff / EP0;
	if (std::abs(4.0 * chi * d) < double(1E-8) * mp.epre[m] * mp.epre[m])
		e[i][j][k] = d / mp.epre[m];
	else if (mp.epre[m] * mp.epre[m] + 4.0 * chi * d >= 0)
		e[i][j][k] = (-mp.epre[m] + sqrt(mp.epre[m] * mp.epre[m] + 4.0 * chi * d)) / (2.0 * chi);
	else
		e[i][j][k] = d / (mp.epre[m] + chi * e[i][j][k]);
}

/**
 * @brief THG: solve E^3 + p E + q = 0, Cardano, p = eps / chi, q = -(dt curl / (eps0 chi) + p E_old + E_old^3)
 * @verbatim
      discriminant D = p^3 / 27 + q^2 / 4
      |D| < 1e-8:   E = -1.5 q / p                          double root; the other root is 3 q / p
      D > 0:        u = cbrt(-q / 2 - sqrt(D)),  E = u - p / (3 u), one real root
      D < 0:        u = 2 sqrt(-p / 3), t = acos(3 q / p / u) / 3, E = u cos(t)
                    first of three roots; the others are u cos(t - 2 pi / 3), u cos(t - 4 pi / 3)
   @endverbatim
 */
void exact_thg(view3<double> e, int m, const double *chi3, double rotdiff, int i, int j, int k)
{
	const auto &mp = fdtd::materials::model_parameters();
	const double chi = chi3[0];
	const double p = mp.epre[m] / chi;
	const double q = -1 * (fdtd::materials::model_spec().dt * rotdiff / (EP0 * chi) + mp.epre[m] / chi * e[i][j][k] + e[i][j][k] * e[i][j][k] * e[i][j][k]);
	const double discr = p * p * p / 27 + q * q / 4;
	if (std::abs(discr) < 1e-8)
		e[i][j][k] = -1.5 * q / p;
	else if (discr > 0)
	{
		const double u = cbrt(-1 * q / 2 - sqrt(discr));
		e[i][j][k] = u - p / (3 * u);
	}
	else
	{
		const double u = 2 * sqrt(-1 * p / 3);
		const double t = acos(3 * q / p / u) / 3;
		e[i][j][k] = u * cos(t);
	}
}
} // namespace

// the diagonal chi element of component a: chi2 xx, yy, zz = [0], [7], [14]; chi3 xxxx, yyyy, zzzz likewise
void fdtd::materials::nonlinear::update_exact_shg(int a, double rotdiff, int m, int i, int j, int k)
{
	exact_shg(fields::fields().E[0][a], m, &chi2_dl[m][7 * a], rotdiff, i, j, k);
}

void fdtd::materials::nonlinear::update_exact_thg(int a, double rotdiff, int m, int i, int j, int k)
{
	exact_thg(fields::fields().E[0][a], m, &chi3_dl[m][7 * a], rotdiff, i, j, k);
}
