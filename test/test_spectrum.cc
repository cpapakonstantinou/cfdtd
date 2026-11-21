// Unit tests of the quadrature, inc/core/quadrature.h, and of the spectral energies of a sampled
// signal, inc/sources/sampled_signal.h, against the closed forms of a Gaussian pulse.
#include "core/quadrature.h"
#include "sources/sampled_signal.h"

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
#define CHECK(c) do { if (!(c)) { std::printf("FAIL %s:%d %s\n", __FILE__, __LINE__, #c); ++failures; } } while (0)

static bool near(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::fabs(b); }

int main()
{
	const double pi = 3.14159265358979323846;

	// panels of n nodes integrate degree 2n - 1 exactly: x^7 on [1, 3] with 3 panels of 4 nodes
	{
		std::vector<double> x;
		std::vector<double> w;
		fdtd::core::gauss_panels(1.0, 3.0, 3, 4, x, w);
		double sum = 0;
		for (std::size_t i = 0; i < x.size(); ++i)
			sum += w[i] * std::pow(x[i], 7);
		CHECK(x.size() == 12 && near(sum, (std::pow(3.0, 8) - 1) / 8, 1e-13));
	}

	// a Gaussian pulse, field FWHM w, carrier f0: exp(-4 ln2 (t - t0)^2 / w^2) cos(2 pi f0 t)
	const double w = 20e-15;
	const double f0 = 2.4e14;
	const double dt = 2e-17;
	const double t0 = 5 * w;
	std::vector<double> s;
	for (int n = 0; n * dt < 2 * t0; ++n)
	{
		const double t = n * dt - t0;
		s.push_back(std::exp(-4 * std::log(2.0) * t * t / (w * w)) * std::cos(2 * pi * f0 * n * dt));
	}
	const fdtd::sources::SampledSignal pulse(s, dt);

	// fluence: w sqrt(pi / (8 ln2)) for a unit envelope, half of it with the carrier
	CHECK(near(pulse.fluence(), 0.5 * w * std::sqrt(pi / (8 * std::log(2.0))), 1e-6));
	// Parseval: the positive frequencies up to Nyquist hold the whole fluence
	CHECK(near(pulse.band_energy(1e9, 0.5 / dt), pulse.fluence(), 1e-6));

	// spectral intensity: a Gaussian of FWHM 2 sqrt(2) ln2 / (pi w), sigma = FWHM / (2 sqrt(2 ln2))
	const double fwhm = 2 * std::sqrt(2.0) * std::log(2.0) / (pi * w);
	const double sigma = fwhm / (2 * std::sqrt(2 * std::log(2.0)));
	CHECK(near(pulse.fwhm(0.5 * f0, 1.5 * f0), fwhm, 1e-3));
	// the share inside f0 +- sigma: erf(1 / sqrt 2)
	CHECK(near(pulse.band_energy(f0 - sigma, f0 + sigma) / pulse.fluence(), std::erf(1 / std::sqrt(2.0)), 1e-4));
	// and the half width that holds it is sigma
	CHECK(near(pulse.half_width_for(f0, std::erf(1 / std::sqrt(2.0)) * pulse.fluence(), 0.5 * f0), sigma, 1e-3));
	// a band beyond the reach of the whole axis
	CHECK(std::isnan(pulse.half_width_for(f0, 2 * pulse.fluence(), 0.5 * f0)));

	// the cube: the third harmonic is sqrt 3 wider, and holds a quarter of the cube's energy at 3 f0
	const fdtd::sources::SampledSignal cube = pulse.power(3);
	CHECK(near(cube.fwhm(2.5 * f0, 3.5 * f0), std::sqrt(3.0) * fwhm, 1e-3));
	const double lobe_1 = cube.band_energy(0.5 * f0, 1.5 * f0);
	const double lobe_3 = cube.band_energy(2.5 * f0, 3.5 * f0);
	CHECK(near(lobe_1 / lobe_3, 9.0, 1e-4)); // cos^3 = (3 cos + cos 3) / 4

	// the rule energy of a fine panel rule equals the band energy
	{
		std::vector<double> f;
		std::vector<double> q;
		fdtd::core::gauss_panels(f0 - sigma, f0 + sigma, 32, 16, f, q);
		CHECK(near(pulse.rule_energy(f, q), pulse.band_energy(f0 - sigma, f0 + sigma), 1e-10));
	}

	if (failures == 0)
		std::printf("test_spectrum: all checks passed\n");
	return failures ? 1 : 0;
}
