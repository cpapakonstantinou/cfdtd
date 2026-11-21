/** \file   src/core/quadrature.cc
 *  \brief  Gauss-Legendre nodes and weights: Newton on P_n from the Chebyshev guesses.
 */
#include "core/quadrature.h"

#include <cmath>

void fdtd::core::gauss_legendre(int n, std::vector<double> &x, std::vector<double> &w)
{
	const double pi = 3.14159265358979323846;
	x.assign(n, 0.0);
	w.assign(n, 0.0);
	for (int i = 0; i < (n + 1) / 2; ++i)
	{
		double z = std::cos(pi * (i + 0.75) / (n + 0.5));
		double dp = 0.0;
		for (int it = 0; it < 100; ++it)
		{
			double p0 = 1.0;
			double p1 = z;
			for (int k = 2; k <= n; ++k)
			{
				const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / k;
				p0 = p1;
				p1 = p2;
			}
			dp = n * (z * p1 - p0) / (z * z - 1.0);
			const double dz = p1 / dp;
			z -= dz;
			if (std::fabs(dz) < 1e-16)
				break;
		}
		x[i] = -z;
		x[n - 1 - i] = z;
		w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * dp * dp);
	}
}

void fdtd::core::gauss_panels(double a, double b, int panels, int n, std::vector<double> &x, std::vector<double> &w)
{
	std::vector<double> t;
	std::vector<double> q;
	gauss_legendre(n, t, q);
	x.clear();
	w.clear();
	const double h = (b - a) / panels;
	for (int p = 0; p < panels; ++p)
	{
		const double mid = a + (p + 0.5) * h;
		for (int i = 0; i < n; ++i)
		{
			x.push_back(mid + 0.5 * h * t[i]);
			w.push_back(0.5 * h * q[i]);
		}
	}
}
