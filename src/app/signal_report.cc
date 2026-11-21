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
 *  Every table is one key: its first line the comma-separated header, then one line per row, so
 *  `grep '^band:'` extracts a CSV table.
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
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <limits>
#include <string>
#include <sstream>
#include <vector>
#include "app/signal_config.h"

namespace
{
void print_usage(std::FILE *to)
{
	std::fprintf(to,
				 "usage: fdtd signal [-s SECTION.KEY=VALUE]... [-t X] [-c X] [-e X] run.ini\n"
				 "Report the source signal, the time steps it needs, its spectrum and the harmonic noise floors.\n"
				 "\n"
				 "  -s, --set SECTION.KEY=VALUE  override a configuration key\n"
				 "  -t, --threshold X            level that ends the signal, relative to its peak; default 1e-3\n"
				 "  -c, --capture X              share of the harmonic the harmonic band must hold; default 0.999\n"
				 "  -e, --efficiency X           expected harmonic efficiency, for the noise margin\n"
				 "  -h, --help                   print this help\n");
}

/// one row of the config table: section, key, value
using ConfigRow = std::array<std::string, 3>;

std::string number(const char *format, double v)
{
	char b[64];
	std::snprintf(b, sizeof b, format, v);
	return b;
}

/// The a priori estimates for the bands of [spectrum]: the share of the pump and of the harmonic the
/// bands hold, the accuracy of the band rule, and the noise efficiency of the harmonic band. The
/// harmonic is that of an instantaneous response, flat over the band: the k-th power of the signal.
void band_report(const fdtd::config::SpectrumSection &s, bool carrier, double f0, double dt, int configured, int span, double capture, double expected,
				 std::vector<ConfigRow> &config_rows)
{
	using fdtd::constants::C0;
	using fdtd::constants::PI2;
	if (s.points < 2)
		return;
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

	auto row = [&](const char *name, int order, double fwhm, double in_band, double rule_error, const char *verdict) {
		const double c = order * f0;
		std::printf("band: %s,%.4f,%.4f,%.4f,%.4f,%s,%s,%s,%s,%s\n", name, s.lambda_min / order * 1e9, s.lambda_max / order * 1e9, order * fa / 1e12,
					order * fb / 1e12, std::isfinite(fwhm) ? number("%.3f", fwhm / 1e12).c_str() : "", std::isfinite(fwhm) ? number("%.2f", C0 * fwhm / (c * c) * 1e9).c_str() : "",
					std::isfinite(in_band) ? number("%.4f", in_band).c_str() : "", std::isfinite(rule_error) ? number("%.2e", rule_error).c_str() : "", verdict);
	};
	std::printf("band: band,lambda_min_nm,lambda_max_nm,f_min_thz,f_max_thz,fwhm_thz,fwhm_nm,in_band,rule_error,verdict\n");
	const double pump_band = pulse.band_energy(fa, fb);
	row("pump", 1, carrier ? pulse.fwhm(0.5 * f0, 1.5 * f0) : pulse.fwhm(0, 2 * fb), pump_band / pulse.fluence(),
		std::fabs(pulse.rule_energy(pump_f, weights) / pump_band - 1), "");
	if (k < 2)
		return;

	const double nan = std::numeric_limits<double>::quiet_NaN();
	if (!carrier)
		row(s.harmonic_band.c_str(), k, nan, nan, nan, "");
	else
	{
		const fdtd::sources::SampledSignal harmonic = pulse.power(k);
		const double lobe = harmonic.band_energy((k - 0.5) * f0, (k + 0.5) * f0);
		const double in_band = harmonic.band_energy(k * fa, k * fb);
		const bool covered = in_band >= capture * lobe;
		row(s.harmonic_band.c_str(), k, harmonic.fwhm((k - 0.5) * f0, (k + 0.5) * f0), in_band / lobe,
			std::fabs(harmonic.rule_energy(harmonic_f, harmonic_w) / in_band - 1), covered ? "covered" : "narrow");
		// a band centred on the carrier in frequency that holds the share asked for
		const double h = covered ? nan : harmonic.half_width_for(k * f0, capture * lobe, 0.5 * f0);
		if (std::isfinite(h))
		{
			config_rows.push_back({"spectrum", "lambda_min", number("%.6e", C0 / (f0 + h / k))});
			config_rows.push_back({"spectrum", "lambda_max", number("%.6e", C0 / (f0 - h / k))});
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
	std::printf("noise: steps,noise_efficiency,rounding_efficiency,expected_efficiency,margin_db\n");
	std::printf("noise: %d,%.3e,%.3e,%s,%s\n", configured, noise, rounding, expected > 0 ? number("%g", expected).c_str() : "",
				expected > 0 ? number("%.1f", 10 * std::log10(expected / std::max(noise, rounding))).c_str() : "");
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
		std::vector<ConfigRow> config_rows;

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
		{
			std::printf("signal: type,dt_s,steps,transit_steps,verdict\n");
			std::printf("signal: %s,%.6e,%d,%d,steady\n", config.source.signal.c_str(), dt, configured, transit);
		}
		else
		{
			e = find_extent(threshold);
			if (e.last_on < 0)
			{
				std::fprintf(stderr, "error: the signal does not end within 50000000 steps\n");
				return 1;
			}
			const int needed = e.last_on + std::max(1, e.last_on / 50) + 1;
			const int recommended = ((needed + 2 * transit) / 10 + 1) * 10;
			std::printf("signal: type,dt_s,steps,transit_steps,peak,peak_step,last_on_step,steps_minimum,steps_recommended,verdict\n");
			std::printf("signal: %s,%.6e,%d,%d,%.6e,%d,%d,%d,%d,%s\n", config.source.signal.c_str(), dt, configured, transit, e.peak, e.peak_step, e.last_on,
						needed, recommended, configured >= needed ? "covered" : "truncated");
			config_rows.push_back({"grid", "time_steps", std::to_string(std::max(configured, recommended))});
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
		std::printf("spectrum: lambda_nm,band,amplitude,relative,flag\n");
		for (std::size_t m = 0; m < omegas.size(); ++m)
		{
			const double rel = top > 0 ? amplitude[m] / top : 0;
			const bool pump = m < pump_points;
			weak += pump && rel < 1e-2;
			std::printf("spectrum: %.4f,%s,%.4e,%.4e,%s\n", 2 * M_PI * fdtd::constants::C0 / omegas[m] * 1e9, pump ? "pump" : "harmonic", amplitude[m], rel,
						!pump ? "raw" : rel < 1e-2 ? "weak" : "ok");
		}

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
		// -250 dB stands for double-precision rounding and below
		auto floor_db = [&](int harmonic, int steps) {
			const double p1 = power_at(f0, steps);
			const double ph = power_at(harmonic * f0, steps);
			return p1 > 0 && ph > 0 ? std::max(-250.0, 10 * std::log10(ph / p1)) : -250.0;
		};
		std::printf("harmonic: fundamental_nm,reference\n");
		std::printf("harmonic: %.4f,%s\n", fdtd::constants::C0 / f0 * 1e9, carrier ? "carrier" : "centre");
		// how the floor moves with the number of steps
		const int complete = steady ? 0 : e.last_on + std::max(1, e.last_on / 50) + 1;
		std::vector<int> trial = {configured / 2, configured, 2 * configured, 4 * configured};
		if (!steady)
			trial.push_back(complete);
		std::sort(trial.begin(), trial.end());
		trial.erase(std::unique(trial.begin(), trial.end()), trial.end());
		std::printf("floor: steps,shg_db,thg_db,note\n");
		for (int n : trial)
			if (n > 0)
				std::printf("floor: %d,%.1f,%.1f,%s\n", n, floor_db(2, n), floor_db(3, n), n == configured ? "configured" : n == complete ? "signal_end" : "");

		band_report(config.spectrum, carrier, f0, dt, configured, steady ? configured : find_extent(1e-8).last_on + 1, capture, expected, config_rows);

		std::printf("config: section,key,value\n");
		for (const ConfigRow &r : config_rows)
			std::printf("config: %s,%s,%s\n", r[0].c_str(), r[1].c_str(), r[2].c_str());
		if (weak)
			std::printf("warning: %d of %zu pump points are below 1e-2 of the strongest, the normalised spectra are noisy there\n", weak, pump_points);
		if (!steady && configured < e.last_on)
			std::printf("warning: the last step cuts the signal off, which raises the harmonic floors\n");
	}
	catch (const std::exception &e)
	{
		std::fprintf(stderr, "error: %s\n", e.what());
		return 1;
	}
	return 0;
}
