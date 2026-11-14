/** \file   boundary/cpml.cc
 *  \brief  CPML grading profiles and psi storage, inc/boundary/cpml.h.
 */
#include "boundary/cpml.h"
#include "core/constants.h"

#include <cmath>
#include <new>

namespace fdtd::boundary
{
namespace
{
/// fills entry m of a profile's table at grading position p (m or m + 1/2) of an N-cell slab
void grade(Profile &pr, cvector<double> (&inv_kappa)[2], int st, int m, double p, int N, double sigma_max, double dt, const Grading &g)
{
	const double u = double(N - p) / double(N);
	const double v = double(std::pow(u, g.exponent));
	const double S = sigma_max * v;
	const double A = g.alpha_max * double(std::pow(double(p) / double(N), g.alpha_exponent));
	const double K = double(1.) + (g.kappa_max - double(1.)) * v;
	const double w = double(std::exp((-1.) * ((S / K + A) * (dt / constants::EP0))));
	pr.b[st][m] = w;
	pr.c[st][m] = g.alpha_max != double(0.) ? S * (w - double(1.)) / (K * (S + A * K)) : (w - double(1.)) / K;
	inv_kappa[st][m] = double(1.) / K;
}
} // namespace

Cpml &cpml()
{
	static Cpml instance;
	return instance;
}

void Cpml::set_profiles(const int cells[3], const double du[3], double dt, const Grading &g)
{
	using namespace index;
	for (int a = 0; a < 3; ++a)
	{
		const int N = cells[a];
		thickness[a] = N;
		// sigma is not divided by the background index, tested 2026-11-06, test/cpml/plane_wave.py:
		// n = 1.5 gains 4.6 dB, n = 3.5 loses 1 dB, a metal background (Re n ~ 0.15) would get 6x sigma
		const double sigma_max = (g.exponent + double(1.)) * double(0.8) / double(377.) / du[a];
		for (int f = 0; f < 2; ++f)
		{
			Profile &pr = profile[a][f];
			auto &ik = grid::stretching().inv_kappa[a][f];
			for (int st = 0; st < 2; ++st)
				ik[st] = vector1<double>(N + 1), pr.b[st] = vector1<double>(N + 1), pr.c[st] = vector1<double>(N + 1);
			for (int m = 0; m <= N; ++m)
				grade(pr, ik, integer, m, double(m), N, sigma_max, dt, g);
			for (int m = 0; m < N; ++m)
				grade(pr, ik, half, m, m + 0.5, N, sigma_max, dt, g);
			ik[half][N] = ik[integer][N];
			pr.b[half][N] = pr.b[integer][N];
			pr.c[half][N] = pr.c[integer][N];
		}
	}
}

bool Cpml::allocate(const grid::Subdomain &s, const bool near[3], const bool far[3], bool harmonic)
{
	this->harmonic = harmonic;
	const int cells[3] = {s.cells_x, s.cells_y, s.cells_z};
	const int inner[3] = {s.inner_dim_x, s.inner_dim_y, s.inner_dim_z};
	try
	{
		for (int a = 0; a < 3; ++a)
		{
			int n[3] = {cells[0], cells[1], cells[2]};
			n[a] -= inner[a];
			for (int face = 0; face < 2; ++face)
			{
				if (!(face == index::near ? near[a] : far[a]))
					continue;
				for (int set = 0; set < (harmonic ? 2 : 1); ++set)
					for (int field = 0; field < 2; ++field)
						for (int c = 0; c < 3; ++c)
							if (c != a)
								psi[set][field][c][a][face] = field3<double>(n[0], n[1], n[2]);
			}
		}
	}
	catch (const std::bad_alloc &)
	{
		return false;
	}
	return true;
}

void Cpml::release()
{
	for (int a = 0; a < 3; ++a)
		for (int f = 0; f < 2; ++f)
			for (int st = 0; st < 2; ++st)
				grid::stretching().inv_kappa[a][f][st].reset(), profile[a][f].b[st].reset(), profile[a][f].c[st].reset();
	for (auto &set : psi)
		for (auto &field : set)
			for (auto &comp : field)
				for (auto &deriv : comp)
					for (auto &p : deriv)
						p.reset();
}
} // namespace fdtd::boundary
