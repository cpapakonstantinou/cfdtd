/** \file   src/sources/sampled_signal.cc
 *  \brief  Spectral energies of a sampled signal: the DFT of the monitors, integrated over bands.
 */
#include "sources/sampled_signal.h"
#include "core/quadrature.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <utility>

namespace
{
const double two_pi = 6.283185307179586;
} // namespace

fdtd::sources::SampledSignal::SampledSignal(std::vector<double> samples, double dt) : _s(std::move(samples)), _dt(dt) {}

double fdtd::sources::SampledSignal::density(double f) const
{
	// the phase factor by recurrence, set exactly every 256 samples so that rounding does not build up
	const std::complex<double> step = std::polar(1.0, two_pi * f * _dt);
	std::complex<double> sum = 0;
	std::complex<double> z = 1;
	for (std::size_t n = 0; n < _s.size(); ++n)
	{
		if (n % 256 == 0)
			z = std::polar(1.0, two_pi * f * _dt * double(n));
		sum += _s[n] * z;
		z *= step;
	}
	return std::norm(sum) * _dt * _dt;
}

double fdtd::sources::SampledSignal::fluence() const
{
	double sum = 0;
	for (double v : _s)
		sum += v * v;
	return sum * _dt;
}

double fdtd::sources::SampledSignal::band_energy(double f0, double f1) const
{
	if (!(f1 > f0))
		return 0;
	// |S|^2 varies on the scale 1 / (N dt): four panels of 16 nodes per such width converge to many digits
	const double scale = 1.0 / (std::max(1, size()) * _dt);
	const int panels = std::clamp(int(std::ceil(4 * (f1 - f0) / scale)), 4, 4096);
	std::vector<double> f;
	std::vector<double> w;
	fdtd::core::gauss_panels(f0, f1, panels, 16, f, w);
	return rule_energy(f, w);
}

double fdtd::sources::SampledSignal::rule_energy(const std::vector<double> &f, const std::vector<double> &w) const
{
	double sum = 0;
	for (std::size_t m = 0; m < f.size(); ++m)
		sum += w[m] * density(f[m]);
	return 2 * sum;
}

double fdtd::sources::SampledSignal::fwhm(double f0, double f1) const
{
	const int n = 4001;
	std::vector<double> f(n);
	std::vector<double> d(n);
	for (int i = 0; i < n; ++i)
	{
		f[i] = f0 + (f1 - f0) * i / (n - 1);
		d[i] = density(f[i]);
	}
	const int top = int(std::max_element(d.begin(), d.end()) - d.begin());
	const double half = 0.5 * d[top];
	auto crossing = [&](int dir) {
		for (int i = top; i + dir >= 0 && i + dir < n; i += dir)
			if (d[i + dir] < half)
				return f[i] + (f[i + dir] - f[i]) * (d[i] - half) / (d[i] - d[i + dir]);
		return std::numeric_limits<double>::quiet_NaN();
	};
	return crossing(1) - crossing(-1);
}

double fdtd::sources::SampledSignal::half_width_for(double fc, double energy, double h_max) const
{
	if (band_energy(fc - h_max, fc + h_max) < energy)
		return std::numeric_limits<double>::quiet_NaN();
	double lo = 0;
	double hi = h_max;
	while (hi - lo > 1e-4 * hi)
	{
		const double mid = 0.5 * (lo + hi);
		(band_energy(fc - mid, fc + mid) < energy ? lo : hi) = mid;
	}
	return hi;
}

fdtd::sources::SampledSignal fdtd::sources::SampledSignal::power(int k) const
{
	std::vector<double> s(_s);
	for (double &v : s)
		v = std::pow(v, k);
	return SampledSignal(std::move(s), _dt);
}
