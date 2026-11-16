/** \file   src/app/signal_report.cc
 *  \brief  Signal settings for an FDTD run: how long the source signal is, how many time steps
 *          the run needs, how well the signal covers the analysed spectrum, and the noise floor a
 *          second or third harmonic has to rise above.
 *
 *    fdtd signal [-s SECTION.KEY=VALUE]... [-t THRESHOLD] run.ini
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
#include "sources/signal.h"
#include "sources/state.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <sstream>
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
				 "\n"
				 "options:\n"
				 "  -s, --set SECTION.KEY=VAL  override a configuration key, repeatable, e.g. -s signal.gaussian.width=20e-15\n"
				 "  -t, --threshold X          the signal counts as on above X of its peak, default 1e-3\n"
				 "  -h, --help                 show this help\n");
}
} // namespace

int fdtd::app::signal_report(int argc, char *argv[])
{
	std::vector<std::string> overrides;
	double threshold = 1e-3;
	static const option long_options[] = {{"set", required_argument, nullptr, 's'}, {"threshold", required_argument, nullptr, 't'},
										  {"help", no_argument, nullptr, 'h'}, {nullptr, 0, nullptr, 0}};
	for (int c; (c = getopt_long(argc, argv, "s:t:h", long_options, nullptr)) != -1;)
		switch (c)
		{
			case 's': overrides.push_back(optarg); break;
			case 't': threshold = std::atof(optarg); break;
			case 'h': print_usage(stdout); return 0;
			default: print_usage(stderr); return 2;
		}
	if (optind + 1 != argc || !(threshold > 0 && threshold < 1))
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
		auto config = fdtd::config::load_sim_config(fdtd::config::build_config(text.str(), file, overrides));
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

		int horizon = configured;
		fdtd::sources::SignalExtent e{};
		if (steady)
			std::printf("signal: steady (no end): the step count is set by the settling time, not by the signal\n");
		else
		{
			// widen the window until the signal ends well inside it
			for (horizon = std::max(configured, 1000);; horizon *= 2)
			{
				e = fdtd::sources::signal_extent(horizon, threshold);
				if (e.peak > 0 && e.last_on < horizon - std::max(1, horizon / 10))
					break;
				if (horizon > 50000000)
				{
					std::printf("signal: no end found within %d steps\n", horizon);
					return 1;
				}
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

		// amplitude of the signal at the analysed frequencies
		const auto omegas = fdtd::config::spectrum_omegas(config.spectrum);
		const int span = steady ? configured : e.last_on + 1;
		std::vector<double> amplitude(omegas.size());
		for (std::size_t m = 0; m < omegas.size(); ++m)
		{
			std::complex<double> sum = 0;
			for (int n = 0; n <= span; ++n)
				sum += fdtd::sources::signal_value(n * dt) * std::exp(std::complex<double>(0, omegas[m] * n * dt));
			amplitude[m] = std::abs(sum) * dt;
		}
		const double top = omegas.empty() ? 0 : *std::max_element(amplitude.begin(), amplitude.end());
		int weak = 0;
		for (std::size_t m = 0; m < omegas.size(); ++m)
		{
			const double rel = top > 0 ? amplitude[m] / top : 0;
			weak += rel < 1e-2;
			std::printf("spectrum: %.4f nm amplitude %.4e relative %.4f%s\n", 2 * M_PI * fdtd::constants::C0 / omegas[m] * 1e9, amplitude[m], rel,
						rel < 1e-2 ? " WEAK" : "");
		}
		if (weak)
			std::printf("spectrum: %d of %zu points below 1e-2 of the strongest: normalised spectra are noisy there (harmonic-band points are expected to be weak)\n",
						weak, omegas.size());

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
					carrier ? "signal carrier" : "[spectrum] centre", config.nonlinear.shg.c_str(), config.nonlinear.thg.c_str());
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
	}
	catch (const std::exception &e)
	{
		std::fprintf(stderr, "error: %s\n", e.what());
		return 1;
	}
	return 0;
}
