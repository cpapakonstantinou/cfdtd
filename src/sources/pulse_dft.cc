/**
 * @file sources/pulse_dft.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-01-19
 * 
 * @brief Initializes the frequencies for monitoring and writes out signal
 * information
 * 
 * @details This file creates the global `analysis_omega` vector which is used for the
 * fourier monitors. It immediately uses this array to output the signal used
 * in the simulation in a text file.
 * 
 * @copyright Copyright 2022
 */
#include "grid/state.h"
#include "grid/subdomain.h"
#include "materials/nonlinear/state.h"
#include "io/text.h"
#include "sources/state.h"
#include "sources/signal.h"
#include "materials/materials.h"
#include "sources/types.h"
#include "core/limits.h"
#include "sources/spec.h"
using fdtd::constants::C0;
using fdtd::constants::INVVACIMP;
using fdtd::constants::I_UNIT;
using fdtd::constants::PI;
using fdtd::constants::PI2;
using fdtd::constants::ROOTPROCESS;
using fdtd::materials::max_slots;
using fdtd::materials::max_slots;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
} // namespace
using namespace fdtd::sources; // transitional: module state, sources/state.h

// File-local state, formerly globals used only in this file.
namespace
{
double input_intensity;                                       ///< Incident intensity |E|^2 at the first analysed frequency in V^2/m^2.
double PhaseNormalization[fdtd::max_frequencies];                ///< Phase of the DFT of the source signal per frequency in rad.
std::complex<double> signal_transform;                                     ///< Scratch: DFT of the source signal at one frequency.
double sum_n;                                                 ///< Scratch sum used to average refractive indices.
int n_max;                                                  ///< Last time step, exclusive, included in the source-signal DFT.
} // namespace

/**
 * @brief Initializes `analysis_omega` and writes out `signal.txt`
 * 
 * @details This function creates `analysis_omega`, a global contruct with all the
 * angular frequencies of interest. The way in which it does this varies
 * depending on the type of simulation being performed. It will also use
 * these values to output the signal used in the simulation in both the 
 * time and frequency domain into reports/source_signal.txt and source_spectrum.txt. 
 */
int fdtd::sources::transform_pulse(int steps)
{
	int n;
	int m;

	// exactly the [spectrum] frequencies, fundamental band, then the harmonic band
	for (m = 0; m < sp.nr_frequencies; ++m)
		fdtd::sources::analysis_omega[m] = spectrum_omega[m];

	const Signal &signal = run_signal();
	using Shape = Signal::Shape;
	if (signal.shape() == Shape::raised_cosine || signal.shape() == Shape::raised_cosine_derivative || signal.shape() == Shape::sinc)
		n_max = int(signal.period() / sp.dt) + 10;
	else if (signal.has_carrier())
		n_max = steps;
	else if (sp.continuous_wave)
	{
		n_max = int(signal.cw_periods() / sp.fundamental_frequency / sp.dt);
		steps = n_max; // a CW run lasts its CW periods
		n_start = int(sp.start_cw / sp.fundamental_frequency / sp.dt);
	}
	else if (sp.cw_pulse)
		n_max = int(signal.cw_periods() / sp.fundamental_frequency / sp.dt) + 10;

	fdtd::array1<std::complex<double>> SignalTransformVec;
	SignalTransformVec = fdtd::vector1<std::complex<double>>(sp.nr_frequencies);
	FILE *sigfile = nullptr;
	if (fdtd::grid::world_rank == ROOTPROCESS)
		sigfile = fdtd::io::open_report(sp.reports_dir + "/source_signal.txt", "signal\tcarrier, one row per step from the DFT start");

	for (m = 0; m < sp.nr_frequencies; ++m)
	{
		signal_transform = std::complex<double>(0.);
		for (int n = n_start; n < n_max; ++n)
		{
			const double tt = n * sp.dt;
			const std::complex<double> w = exp(I_UNIT * fdtd::sources::analysis_omega[m] * tt);
			signal_transform += fdtd::sources::signal_value(tt) * w;
			if (sigfile && m == 0)
				fprintf(sigfile, "%e\t%e\n", fdtd::sources::signal_value(tt), cos(2.0 * PI * signal.carrier_frequency() * tt));
		}
		fdtd::sources::amplitude_normalization[m] = double(sqrt(double(signal_transform.real() * signal_transform.real() + signal_transform.imag() * signal_transform.imag())));
		PhaseNormalization[m] = double(atan2(signal_transform.imag(), signal_transform.real()));
		fdtd::sources::dft_normalization[m] = std::complex<double>(signal_transform.real(), signal_transform.imag());
		SignalTransformVec[m] = signal_transform;

	}
	if (sigfile)
		fclose(sigfile);

	if (fdtd::grid::world_rank == ROOTPROCESS && !sp.grating)
	{
		fdtd::io::log("sources", "Time steps required for excitation: %d, DFT calculated from step %d\n", n_max, n_start);

		if (FILE *string_DFT = fdtd::io::open_report(sp.reports_dir + "/source_spectrum.txt", "LAMBDA[nm]\t DFT_module\t DFT_phase\t N_samples\t DFT_re\t DFT_im"))
		{
			for (m = 0; m < sp.nr_frequencies; ++m)
			{
				fprintf(string_DFT, "%e\t%e\t%e\t%d\t%e\t%e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), fdtd::sources::amplitude_normalization[m], PhaseNormalization[m], n_max - n_start, fdtd::sources::dft_normalization[m].real(), fdtd::sources::dft_normalization[m].imag());
			}

			fclose(string_DFT);
		}
	}

	if (sp.nr_frequencies == 2 && (sp.shg || sp.thg))
	{
		input_intensity = fdtd::sources::amplitude_normalization[0] * fdtd::sources::amplitude_normalization[0];					  // V^2/m^2
		input_power = input_intensity * INVVACIMP * fdtd::materials::refractive_index(sp.medium[0], fdtd::sources::analysis_omega[0]).real(); // W/m^2
		;
		fdtd::sources::amplitude_normalization[1] = double(1.); // amplitude_normalization[0];
		fdtd::sources::dft_normalization[1] = std::complex<double>(1, 0);	  // dft_normalization[0];
		PhaseNormalization[1] = double(0.);	  // PhaseNormalization[0];
	}

	if ((sp.thg || sp.shg || (sp.hydrodynamics && !sp.harmonic)) && (sp.nr_frequencies > 2))
	{
		int max_index = 0.5 * sp.nr_frequencies;
		input_intensity = double(0.);
		input_power = double(0.);
		for (m = 0; m < max_index; ++m)
		{
			input_intensity += 2 * sp.dt * sp.dt * fdtd::sources::amplitude_normalization[m] * fdtd::sources::amplitude_normalization[m] * fundamental_step_frequency;															// s*V^2/m^2
			input_power += 2 * sp.dt * sp.dt * fdtd::sources::amplitude_normalization[m] * fdtd::sources::amplitude_normalization[m] * fundamental_step_frequency * INVVACIMP * fdtd::materials::refractive_index(sp.medium[0], fdtd::sources::analysis_omega[m]).real(); // J/m^2
			fdtd::sources::amplitude_normalization[m + max_index] = double(1.);																											// amplitude_normalization[m];
			fdtd::sources::dft_normalization[m + max_index] = std::complex<double>(1, 0);																												// dft_normalization[m];
			PhaseNormalization[m + max_index] = double(0.);																												// PhaseNormalization[m];
		}
		;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for (n = 0; n < fdtd::materials::slot_count(); ++n)
	{
		sum_n = double(0.);
		for (m = 0; m < sp.nr_frequencies; ++m)
			sum_n += fdtd::materials::refractive_index(n, fdtd::sources::analysis_omega[m]).real();
		if (fdtd::grid::world_rank == ROOTPROCESS && sum_n > 0)
			fdtd::io::log("n_avg", "material %d %s index %f", n, fdtd::materials::slot_name(n).c_str(), sum_n / sp.nr_frequencies);
		if (fdtd::grid::world_rank == ROOTPROCESS && (sp.shg || sp.thg) && fdtd::materials::nonlinear::nl[n])
			fdtd::io::log("sources", "material %d %s is nonlinear\n", n, fdtd::materials::slot_name(n).c_str());
	}

	return steps;
}
