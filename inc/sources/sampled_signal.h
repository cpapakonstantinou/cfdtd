/**
 * @file sources/sampled_signal.h
 * @brief Spectral energies of a sampled signal, the a priori estimates of `fdtd signal`
 *
 * The spectrum is the DFT as the monitors take it, S(f) = dt sum_n s_n exp(2 pi i f n dt), and the
 * energies are those of the harmonic efficiency: 2 times the integral of |S(f)|^2 over a band of
 * positive frequencies, so that the whole positive axis gives the fluence sum_n s_n^2 dt.
 */
#ifndef FDTD_SOURCES_SAMPLED_SIGNAL_H
#define FDTD_SOURCES_SAMPLED_SIGNAL_H

#include <vector>

namespace fdtd::sources
{
class SampledSignal
{
public:
	/// samples s_n = s(n dt), n = 0 .. N - 1, of a real signal
	SampledSignal(std::vector<double> samples, double dt);

	double dt() const { return _dt; }
	int size() const { return int(_s.size()); }

	/// |S(f)|^2
	double density(double f) const;
	/// sum_n s_n^2 dt
	double fluence() const;
	/// 2 times the integral of |S(f)|^2 over [f0, f1], Gauss-Legendre panels fine enough for the signal's length
	double band_energy(double f0, double f1) const;
	/// 2 sum_m w_m |S(f_m)|^2, a rule with nodes f in Hz and weights w in Hz
	double rule_energy(const std::vector<double> &f, const std::vector<double> &w) const;
	/// full width at half maximum of |S(f)|^2 around its largest value in [f0, f1]; NaN when a side stays above half inside
	double fwhm(double f0, double f1) const;
	/// the smallest h for which band_energy(fc - h, fc + h) reaches `energy`, searched up to h_max; NaN when even h_max does not
	double half_width_for(double fc, double energy, double h_max) const;
	/// the samples raised to the power k
	SampledSignal power(int k) const;

private:
	std::vector<double> _s;
	double _dt;
};
} // namespace fdtd::sources

#endif // FDTD_SOURCES_SAMPLED_SIGNAL_H
