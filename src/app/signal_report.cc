/** \file   src/app/signal_report.cc
 *  \brief  Signal settings for an FDTD run: how long the source signal is, how many time steps
 *          the run needs, how well the signal covers the analysed spectrum, and the noise floor a
 *          second or third harmonic has to rise above; for a harmonic band, how much of the pump
 *          and of the harmonic the bands hold and the noise in the harmonic band.
 *
 *    fdtd signal [-s SECTION.KEY=VALUE]... [-t THRESHOLD] [-c CAPTURE] [-e EFFICIENCY] run.ini
 *
 *  It builds the signal exactly as the solver does, the `sources` module, from the same
 *  configuration, and samples it. No MPI, no fields: it runs on a login node in well under a second.
 *  Every output line starts with a key, so `grep '^signal:'` or `grep '^config:'` extracts a part.
 **/
#include "app/signal_report.h"
#include "app/run_config.h"
#include "app/run_state.h"
#include "app/setup.h"
#include "config/apply.h"
#include "core/constants.h"
#include "config/sim_config.h"
#include "sources/sampled_signal.h"
#include "sources/signal.h"
#include "sources/state.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <limits>
#include <sstream>
#include <vector>
#include "app/signal_config.h"

namespace
{
void print_usage(std::FILE *to)
{
	std::fprintf(to,
				 "usage: fdtd signal [options] run.ini\n"
				 "\n"
				 "Reports the length of the source signal of a run, the time steps it needs and the\n"
				 "signal amplitude at the analysed wavelengths, and the pump's own power at the second and\n"
				 "third harmonic, the floor a harmonic signal must exceed, for the configured time steps.\n"
				 "With a [spectrum] harmonic_band, the share of the pump and of the harmonic inside the bands,\n"
				 "the accuracy of the band rule and the noise efficiency of the harmonic band.\n"
				 "\n"
				 "options:\n"
				 "  -s, --set SECTION.KEY=VAL  override a configuration key, repeatable, e.g. -s signal.gaussian.width=20e-15\n"
				 "  -t, --threshold X          the signal counts as on above X of its peak, default 1e-3\n"
				 "  -c, --capture X            share of the harmonic the band should hold, default 0.999\n"
				 "  -e, --efficiency X         expected harmonic efficiency, for the margin above the noise\n"
				 "  -h, --help                 show this help\n");
}

/// The a priori estimates for the bands of [spectrum]: the share of the pump and of the harmonic the
/// bands hold, the accuracy of the band rule, and the noise efficiency of the harmonic band. The
/// harmonic is that of an instantaneous response, flat over the band: the k-th power of the signal.
void band_report(const fdtd::config::SpectrumSection &s, bool carrier, double f0, double dt, int configured, int span, double capture, double expected)
{
	using fdtd::constants::C0;
	using fdtd::constants::PI2;
	if (s.points < 2)
	{
		std::printf("band: a single wavelength, no band\n");
		return;
	}
	const int k = fdtd::config::harmonic_order(s.harmonic_band);
	const double fa = C0 / s.lambda_max;
	const double fb = C0 / s.lambda_min;
	const std::vector<double> omegas = fdtd::config::spectrum_omegas(s, C0, PI2);
	const std::vector<double> weights = fdtd::config::spectrum_weights(s, C0, PI2);
	std::vector<double> pump_f;
	std::vector<double> harmonic_f;
	std::vector<double> harmonic_w;
	for (int m = 0; m < s.points; ++m)
	{
		pump_f.push_back(omegas[m] / PI2);
		harmonic_f.push_back(k * omegas[m] / PI2);
		harmonic_w.push_back(k * weights[m]);
	}

	// the whole signal, every stride-th step: 32 samples per period of the highest frequency analysed
	const double f_top = std::max(k * fb, carrier ? (k + 0.5) * f0 : 2 * fb);
	const int stride = std::max(1, int(1 / (32 * f_top * dt)));
	std::vector<double> samples;
	for (int n = 0; n < span; n += stride)
		samples.push_back(fdtd::sources::signal_value(n * dt));
	const fdtd::sources::SampledSignal pulse(std::move(samples), stride * dt);

	std::printf("band: pump %.4f to %.4f nm, %.4f to %.4f THz, %d points, %s spacing\n", s.lambda_min * 1e9, s.lambda_max * 1e9, fa / 1e12, fb / 1e12,
				s.points, s.spacing.c_str());
	const double pump_fwhm = carrier ? pulse.fwhm(0.5 * f0, 1.5 * f0) : pulse.fwhm(0, 2 * fb);
	if (std::isfinite(pump_fwhm))
		std::printf("band: pump_fwhm %.3f THz, %.2f nm at the fundamental (spectral intensity of the signal)\n", pump_fwhm / 1e12, C0 * pump_fwhm / (f0 * f0) * 1e9);
	else
		std::printf("band: pump_fwhm wider than the search range\n");
	const double pump_band = pulse.band_energy(fa, fb);
	std::printf("band: pump_in_band %.4f (share of the input fluence; the efficiency takes its input from the time signal, so this limits only the linear spectra)\n",
				pump_band / pulse.fluence());
	std::printf("band: pump_rule_error %.2e (the %d-point rule of the band against a converged integral)\n",
				std::fabs(pulse.rule_energy(pump_f, weights) / pump_band - 1), s.points);
	if (k < 2)
		return;

	std::printf("band: harmonic %s %.4f to %.4f nm, %.4f to %.4f THz\n", s.harmonic_band.c_str(), s.lambda_min / k * 1e9, s.lambda_max / k * 1e9, k * fa / 1e12,
				k * fb / 1e12);
	if (!carrier)
		std::printf("band: harmonic share not estimated: the signal has no carrier, its harmonics overlap\n");
	else
	{
		const fdtd::sources::SampledSignal harmonic = pulse.power(k);
		const double lobe = harmonic.band_energy((k - 0.5) * f0, (k + 0.5) * f0);
		const double in_band = harmonic.band_energy(k * fa, k * fb);
		const double harmonic_fwhm = harmonic.fwhm((k - 0.5) * f0, (k + 0.5) * f0);
		if (std::isfinite(harmonic_fwhm))
			std::printf("band: harmonic_fwhm %.3f THz, %.2f nm at %d f (spectral intensity of the signal to the power %d)\n", harmonic_fwhm / 1e12,
						C0 * harmonic_fwhm / (k * f0 * k * f0) * 1e9, k, k);
		std::printf("band: harmonic_in_band %.4f (share of the harmonic energy near %d f; the harmonic efficiency is low by this factor)\n", in_band / lobe, k);
		std::printf("band: harmonic_rule_error %.2e (the %d-point rule of the harmonic band against a converged integral)\n",
					std::fabs(harmonic.rule_energy(harmonic_f, harmonic_w) / in_band - 1), s.points);
		if (in_band >= capture * lobe)
			std::printf("band: verdict covered (the band holds at least %g of the harmonic)\n", capture);
		else
		{
			std::printf("band: verdict NARROW: %.1f %% of the harmonic falls outside the band; widen the band or lengthen the pulse\n", 100 * (1 - in_band / lobe));
			const double h = harmonic.half_width_for(k * f0, capture * lobe, 0.5 * f0);
			if (std::isfinite(h))
			{
				std::printf("band: the band below is centred on the carrier in frequency and holds %g of the harmonic\n", capture);
				std::printf("config: [spectrum] lambda_min = %.4e\n", C0 / (f0 + h / k));
				std::printf("config: [spectrum] lambda_max = %.4e\n", C0 / (f0 - h / k));
			}
		}
	}

	// the pump's own energy in the harmonic band: the DFT over the configured steps, at the solver's nodes and weights
	std::vector<double> run(std::max(configured, 0));
	for (int n = 0; n < configured; ++n)
		run[n] = fdtd::sources::signal_value(n * dt);
	const fdtd::sources::SampledSignal window(std::move(run), dt);
	const double fluence = window.fluence();
	const double noise = fluence > 0 ? window.rule_energy(harmonic_f, harmonic_w) / fluence : 0;
	// rounding of N terms of the DFT sum, a random walk of eps |s_n|: energy 2 eps^2 dt^2 sum s_n^2 per Hz
	const double eps = 0.5 * std::numeric_limits<double>::epsilon();
	const double rounding = 2 * eps * eps * dt * k * (fb - fa);
	std::printf("band: noise_efficiency %.3e (the pump's own energy in the harmonic band over the input fluence, %d steps; for a linear response of 1 a harmonic efficiency near this is the pump)\n",
				noise, configured);
	std::printf("band: rounding_efficiency %.3e (double rounding of the DFT sums alone, %d steps)\n", rounding, configured);
	if (expected > 0)
		std::printf("band: margin %.1f dB (expected efficiency %g over the larger of the two)\n", 10 * std::log10(expected / std::max(noise, rounding)), expected);
}
} // namespace

int fdtd::app::signal_report(int argc, char *argv[])
{
	std::vector<std::string> overrides;
	double threshold = 1e-3;
	double capture = 0.999;
	double expected = 0;
	static const option long_options[] = {{"set", required_argument, nullptr, 's'}, {"threshold", required_argument, nullptr, 't'},
										  {"capture", required_argument, nullptr, 'c'}, {"efficiency", required_argument, nullptr, 'e'},
										  {"help", no_argument, nullptr, 'h'}, {nullptr, 0, nullptr, 0}};
	for (int c; (c = getopt_long(argc, argv, "s:t:c:e:h", long_options, nullptr)) != -1;)
		switch (c)
		{
			case 's': overrides.push_back(optarg); break;
			case 't': threshold = std::atof(optarg); break;
			case 'c': capture = std::atof(optarg); break;
			case 'e': expected = std::atof(optarg); break;
			case 'h': print_usage(stdout); return 0;
			default: print_usage(stderr); return 2;
		}
	if (optind + 1 != argc || !(threshold > 0 && threshold < 1) || !(capture > 0 && capture < 1) || expected < 0)
	{
		print_usage(stderr);
		return 2;
	}
	const std::string file = argv[optind];
	std::ifstream in(file);
	if (!in)
	{
		std::fprintf(stderr, "error: cannot read configuration file '%s'\n", file.c_str());
		return 1;
	}
	std::stringstream text;
	text << in.rdbuf();

	try
	{
		overrides.push_back("source.allow_truncated_signal=true"); // this tool reports instead of failing
		auto config = fdtd::config::load_sim_config(fdtd::config::build_config(fdtd::config::with_catalog(text.str(), file, overrides), file, overrides));
		fdtd::config::validate_sim_config(config);
		fdtd::config::apply_sim_config(config, 1, true);
		fdtd::app::set_time_step();
		fdtd::app::configure_signal();

		const auto &rc = fdtd::app::run_config();
		const double dt = rc.time.dt;
		const int configured = fdtd::app::run_state().nr_max_time_steps;
		const bool steady = fdtd::sources::run_signal().steady();

		// time for light to cross the longest side of the domain in vacuum
		const double side = std::max({rc.grid.nr_cells_tot_x * rc.grid.dx, rc.grid.nr_cells_tot_y * rc.grid.dy, rc.grid.nr_cells_tot_z * rc.grid.dz});
		const int transit = int(std::ceil(side / (fdtd::constants::C0 * dt)));

		std::printf("signal: type %s\n", config.source.signal.c_str());
		std::printf("signal: dt %.6e s\n", dt);
		std::printf("signal: steps_configured %d (%.4e s)\n", configured, configured * dt);
		std::printf("signal: transit_steps %d (vacuum, longest side %.4e m)\n", transit, side);

		// widen the window until the signal ends well inside it; last_on -1 when no end is found
		auto find_extent = [&](double level) {
			for (int horizon = std::max(configured, 1000); horizon <= 50000000; horizon *= 2)
			{
				const auto found = fdtd::sources::signal_extent(horizon, level);
				if (found.peak > 0 && found.last_on < horizon - std::max(1, horizon / 10))
					return found;
			}
			return fdtd::sources::SignalExtent{0, 0, -1};
		};
		fdtd::sources::SignalExtent e{};
		if (steady)
			std::printf("signal: steady (no end): the step count is set by the settling time, not by the signal\n");
		else
		{
			e = find_extent(threshold);
			if (e.last_on < 0)
			{
				std::printf("signal: no end found within 50000000 steps\n");
				return 1;
			}
			const int needed = e.last_on + std::max(1, e.last_on / 50) + 1;
			const int recommended = ((needed + 2 * transit) / 10 + 1) * 10;
			std::printf("signal: peak %.6e at step %d (%.4e s)\n", e.peak, e.peak_step, e.peak_step * dt);
			std::printf("signal: last_step_on %d (%.4e s, above %g of the peak)\n", e.last_on, e.last_on * dt, threshold);
			std::printf("signal: steps_minimum %d (the solver's guard: the signal ends before the last 2 %% of the run)\n", needed);
			std::printf("signal: steps_recommended %d (signal + 2 transits; resonant structures ring longer)\n", recommended);
			std::printf("signal: verdict %s\n", configured >= needed ? "covered" : "TRUNCATED: the run would stop at the guard");
			std::printf("config: [grid] time_steps = %d\n", std::max(configured, recommended));
		}

		// amplitude of the signal at the analysed frequencies, the DFT over the run's steps as the monitors take it
		const auto omegas = fdtd::config::spectrum_omegas(config.spectrum, fdtd::constants::C0, fdtd::constants::PI2);
		const std::size_t pump_points = std::size_t(config.spectrum.points);
		std::vector<double> amplitude(omegas.size());
		for (std::size_t m = 0; m < omegas.size(); ++m)
		{
			std::complex<double> sum = 0;
			for (int n = 0; n < configured; ++n)
				sum += fdtd::sources::signal_value(n * dt) * std::exp(std::complex<double>(0, omegas[m] * n * dt));
			amplitude[m] = std::abs(sum) * dt;
		}
		// the linear spectra divide by the pump points; the harmonic band is the raw DFT, not normalised
		const double top = omegas.empty() ? 0 : *std::max_element(amplitude.begin(), amplitude.begin() + std::min(pump_points, omegas.size()));
		int weak = 0;
		for (std::size_t m = 0; m < omegas.size(); ++m)
		{
			const double rel = top > 0 ? amplitude[m] / top : 0;
			const bool pump = m < pump_points;
			weak += pump && rel < 1e-2;
			std::printf("spectrum: %.4f nm amplitude %.4e relative %.4e%s\n", 2 * M_PI * fdtd::constants::C0 / omegas[m] * 1e9, amplitude[m], rel,
						!pump ? " harmonic (raw DFT, not normalised)" : rel < 1e-2 ? " WEAK" : "");
		}
		if (weak)
			std::printf("spectrum: %d of %zu pump points below 1e-2 of the strongest: the normalised spectra are noisy there\n", weak, pump_points);

		// The pump's own spectral power at 2 f and 3 f, from a DFT over the run's steps as the solver's
		// monitors take it: the floor below which a harmonic is not seen. Cutting the signal off, or,
		// for a steady signal, a short run, widens its spectrum and raises this floor.
		const bool carrier = fdtd::sources::run_signal().has_carrier();
		const double f0 = carrier ? fdtd::constants::C0 / rc.source.center_freq : rc.source.fundamental_frequency;
		auto power_at = [&](double f, int steps) {
			std::complex<double> sum = 0;
			const double w = 2 * M_PI * f * dt;
			for (int n = 0; n < steps; ++n)
				sum += fdtd::sources::signal_value(n * dt) * std::exp(std::complex<double>(0, w * n));
			return std::norm(sum);
		};
		auto floor_db = [&](int harmonic, int steps) {
			const double p1 = power_at(f0, steps);
			const double ph = power_at(harmonic * f0, steps);
			return p1 > 0 && ph > 0 ? 10 * std::log10(ph / p1) : -INFINITY;
		};
		auto show = [](double db) {
			char b[64];
			if (db < -250)
				std::snprintf(b, sizeof b, "< -250 dB (double-precision rounding)");
			else
				std::snprintf(b, sizeof b, "%.1f dB", db);
			return std::string(b);
		};
		std::printf("harmonic: fundamental %.4f nm (%s); [nonlinear] shg = %s, thg = %s\n", fdtd::constants::C0 / f0 * 1e9,
					carrier ? "signal carrier" : "[spectrum] centre", config.nonlinear.shg ? "true" : "false", config.nonlinear.thg ? "true" : "false");
		for (int h = 2; h <= 3; ++h)
			std::printf("harmonic: %s %.4f nm floor %s (pump power at %d f relative to f, %d steps)\n", h == 2 ? "shg" : "thg",
						fdtd::constants::C0 / (h * f0) * 1e9, show(floor_db(h, configured)).c_str(), h, configured);
		// how the floor moves with the number of steps
		const int complete = steady ? 0 : e.last_on + std::max(1, e.last_on / 50) + 1;
		std::vector<int> trial = {configured / 2, configured, 2 * configured, 4 * configured};
		if (!steady)
			trial.push_back(complete);
		std::sort(trial.begin(), trial.end());
		trial.erase(std::unique(trial.begin(), trial.end()), trial.end());
		for (int n : trial)
			if (n > 0)
				std::printf("harmonic: steps %d shg_floor %s thg_floor %s%s\n", n, show(floor_db(2, n)).c_str(), show(floor_db(3, n)).c_str(),
							n == configured ? " (configured)" : n == complete ? " (signal above the threshold covered)" : "");
		if (!steady && configured < e.last_on)
			std::printf("harmonic: the signal is cut off at the last step: its spectrum spreads and the floors are raised; run to the signal's end\n");

		band_report(config.spectrum, carrier, f0, dt, configured, steady ? configured : find_extent(1e-8).last_on + 1, capture, expected);
	}
	catch (const std::exception &e)
	{
		std::fprintf(stderr, "error: %s\n", e.what());
		return 1;
	}
	return 0;
}
