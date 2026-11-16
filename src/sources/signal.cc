/**
 * @file sources/signal.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-01-19
 * 
 * @brief houses the envelopes, carrier and signal functions.
 * 
 * @copyright Copyright 2022
 */
#include <algorithm>
#include <cmath>
#include <string>
#include "sources/state.h"
#include "sources/signal.h"
#include "sources/kinds.h"
using fdtd::constants::C0;
using fdtd::constants::I_UNIT;
using fdtd::constants::PI;
using fdtd::constants::PI2;

namespace
{
} // namespace
using namespace fdtd::sources; // transitional: module state, sources/state.h

namespace
{
fdtd::sources::Signal the_signal; ///< the run's signal (configure_signal)
} // namespace

const fdtd::sources::Signal &fdtd::sources::run_signal() { return the_signal; }
void fdtd::sources::set_run_signal(const Signal &s) { the_signal = s; }

/**
 * @brief All available envelopes for the signal
 * 
 * @param type [in] integer value corresponding to envelope type
 * @param t [in] time envelope is being evaluated at
 * @param parameter1 [in] first input parameter for signal
 * @param parameter2 [in] second input parameter for signal
 * @param ePeak [in] peak amplitude of signal
 * @return double 
 */
double envelopes(int type, const double &t, const double &parameter1, const double &parameter2, const double &ePeak )
{
	double x;
	double y;
	double tt;
	double period;
	double t_0;
	double delta_t;

	switch (type) {
		// Raised Cosine Source
		case envelope::raised_cosine:
			if ((t > double(0.)) && (t < parameter2))
			{
				tt = double(parameter1 * t);
				y = 1. - cos(tt);
				return ePeak * double(y * y * y / 8.);
			}
			else
				return double(0.);
		// Derived Raised cosine, AKA Derivative of Raised Cosine
		case envelope::raised_cosine_derivative:
			if ((t > double(0.)) && (t < parameter2))
			{
				tt = double(parameter1 * t);
				x = sin(tt);
				y = 1. - cos(tt);
				return ePeak * double(x * y * y / 4.);
			}
			else
				return double(0.);
		// Rectangular Envelope, Top Hat
		case envelope::rectangular:
			if ((t > double(0.) + 100 * parameter1) && (t < parameter2 + 100 * parameter1))
				return ePeak * double(1.);
			else
				return double(0.);
		// CW, Continuous Wave, -- A Very LONG envelope
		case envelope::gaussian_periods: // parameter1: the fundamental frequency
			period = double(1 / parameter1);
			t_0 = double(0.5) * parameter2 * period;
			delta_t = double(0.1) * parameter2 * period;

			if ((t > double(0.)) && (t < parameter2 * period))
			{
				return ePeak * exp(-(t - t_0) * (t - t_0) / (delta_t * delta_t));
			}
			else
				return double(0.);
		// Long Sine Wave
		case envelope::arctan_ramp:
			period = double(1 / parameter1);
			t_0 = 0.5 * parameter2 * period;
			return ePeak * (1 / PI) * (0.5 * PI + atan((t - t_0) / period));
		// Gaussian Envelope
		case envelope::gaussian:
			return ePeak * exp(-2.7726 * (double(t - parameter1) / double(parameter2)) * (double(t - parameter1) / double(parameter2)));
		// Undefined
		case envelope::linear_ramp:
			if (t < parameter2)
			{
				return ePeak * t / parameter2 * parameter1;
			}
			else
				return ePeak;
		// No envelope case
		default:
			return double(1.);
	}

}

/**
 * @brief All available carriers for the signal
 * 
 * @param type [in] integer value corresponding to carrier type
 * @param t [in] time the carrier is being evaluated at
 * @param wavelength [in] wavelength of the carrier
 * @param offset [in] time/phase delay for the carrier
 * @return double 
 */
double carriers(int type, const double &t, const double &wavelength, const double &offset)
{
	double carrier;
	switch (type)
	{
		// Cosine Carrier
		case carrier::cosine:
			carrier = cos(2.0 * PI * wavelength * (t - offset));
			return carrier;
		// Sine Carrier
		case carrier::sine:
			carrier = sin(wavelength * (t - offset));
			return carrier;
		// No Carrier
		default:
			return (double(1));
	}
}

std::complex<double> complex_carriers(int type, const double& t, const double& wavelength, const double& chirp1, const double& chirp2, const std::complex<double>& iUnit, const double& delay, const double& width) {
	std::complex<double> carrier;
	switch (type)
	{
		// Undefined
		case carrier::chirped:
			carrier = std::complex<double>(exp(-iUnit * 2.0 * PI * wavelength * t) *
				exp(iUnit * chirp1 * 2.7726 * (double(t - delay) / double(width)) * (double(t - delay) / double(width))) *
				exp(iUnit * chirp2 * 2.7726 * (double(t - delay) / double(width)) * (double(t - delay) / double(width)) * (double(t - delay) / double(width))));
			return carrier;
		// No Carrier
		default:
			return (std::complex<double>(1));
	}
}

/**
 * @brief Signal function, Electric Field Source
 * 
 * @param t time in seconds
 * @return double signal value in V/m
 * 
 * @details This function uses global vars to determine signal at a given time
 * 't'. This function handles: \n
 * - Raised Cosine \n
 * - Derived Raised Cosine \n
 * - CW Pulse \n
 * - Continuous Wave \n
 * - Modulated Gaussian \n
 * - Constant field \n
 * 
 */
double fdtd::sources::Signal::value(double t) const
{
	double period;
	double k_period;

	// Raised Cosine, no carrier
	if (_shape == Shape::raised_cosine)
		return (envelopes(envelope::raised_cosine,t,_k, _period,_e_peak));

	// Derived Raised Cosine, no carrier
	else if (_shape == Shape::raised_cosine_derivative)
		return (envelopes(envelope::raised_cosine_derivative, t, _k, _period, _e_peak));

	// Rectangular fdtd::sources::signal_value(top hat, no carrier)
	else if (_shape == Shape::rectangular)
		return (envelopes(envelope::rectangular, t, _dt, _period, _e_peak));

	// Continuous Wave pulse, for SHG an THG
	else if (_shape == Shape::cw_pulse)
	{
		period = double(1 / _fundamental);
		k_period = double(PI2 / period);
		return (envelopes(envelope::gaussian_periods, t, _fundamental, _cw_periods, _e_peak) * carriers(carrier::sine, t, k_period, double(double(0.5) * _cw_periods * period)));
	}

	// Long Sine Wave
	else if (_shape == Shape::continuous_wave)
	{
		k_period = double(PI2 * _fundamental);
		return (envelopes(envelope::arctan_ramp, t, _fundamental, _start_cw, _e_peak) * carriers(carrier::sine, t, k_period, double(0)));
	}

	// Modulated Gaussian, can be chirped
	else if (_shape == Shape::gaussian)
	{
			if(_chirp == 0 && _chirp_2 == 0)
			{
				double pulse;
				pulse = envelopes(envelope::gaussian, t, _delay, _width, _e_peak)*carriers(carrier::cosine,t,_f_0,double(0));
				return pulse;
			}
			else
			{
				std::complex<double> pulse;
				pulse = std::complex<double>( envelopes(envelope::gaussian, t, _delay, _width, _e_peak) * complex_carriers(carrier::chirped, t, _f_0,_chirp,_chirp_2,I_UNIT,_delay,_width )) ;
				return pulse.real();
			}

	}

	// Modulated Gaussian 2, cannot be chirped
	else if (_shape == Shape::gaussian_alt)
	{
		double pulse;
		double pulse2;
		if (t < 3 * (_delay + 2 * _width))
		{
			pulse = envelopes(envelope::gaussian, t, _delay, _width, _e_peak) * carriers(carrier::cosine, t, _f_0, double(0));
			pulse2 = envelopes(envelope::gaussian, t, _delay, _width, _e_peak) * carriers(carrier::cosine, t, _f_02, double(0));
		}
		else
		{
			pulse = 0.0;
			pulse2 = 0.0;
		}
		return pulse + pulse2;
	}

	// Contant Field, build up and hold
	else if (_shape == Shape::constant)
	{
		double pulse;
		if (_shape == Shape::constant)
		{
			pulse = envelopes(envelope::linear_ramp,t,_dt, 1000.0 * _dt,_e_peak);
		}
		return pulse;
	}

	else if (_shape == Shape::sinc)
	{	
		/* 
		See Emadeldeen Hassan and Antonio Calà Lesina, 
		"Topology optimization of dispersive plasmonic nanostructures in the time-domain," 
		Opt. Express 30, 19557-19572, 2022
		*/
		if ((t > double(0.)) && (t < _period))
		{	
			double lc = 520e-9;  // 
			double fc = C0 / lc; // Central frquency
			double L = 0.5 * _period;
			int lobes = 7; // broadband 350 - 1000, hald-bandwith bounds
			double carrierSignal = cos(2 * PI * fc * (t - L));
			double truncatedSinc = sin(lobes * PI * (t - L) / L) / (PI * (t - L));
			double hanningWindow = pow(cos(PI * (t - L) / (2 * L)), 2);	
			double normConstant = lobes / L; // #L'Hospital <3

			return carrierSignal * truncatedSinc * hanningWindow / normConstant;
		}
		else
			return double(0.);
	}

	// No source -- return 0
	else
		return double(0.);
}

/**
 * @brief the analytical derivatives of all the envelopes functions.
 * 
 * @param type [in] integer value corresponding to envelope type
 * @param t [in] time envelope is being evaluated at
 * @param parameter1 [in] first input parameter for signal
 * @param parameter2 [in] second input parameter for signal
 * @param ePeak [in] peak amplitude of signal
 * @return double 
 */
double d_envelopes(int type, const double &t, const double &parameter1, const double &parameter2,const double &ePeak )
{
	double x;
	double y;
	double tt;
	double period;
	double t_0;
	double delta_t;
	double exponent;

	switch (type) {
		// Raised Cosine Analytical Derivative
		case envelope::raised_cosine:
			if ((t > double(0.)) && (t < parameter2))
			{
				tt = double(parameter1 * t);
				y = 1. - cos(tt);
				return double(3. * ePeak * parameter1 * sin(tt) * y * y  / 8.);
			}
			else
			{
				return double(0.);
			}

		// Derived Raised Cosine Analytical Derivative
		case envelope::raised_cosine_derivative:
			if ((t > double(0.)) && (t < parameter2))
			{
				tt = double(parameter1 * t);
				x = sin(tt);
				y = 1. - cos(tt);
				return double(ePeak * parameter1 * (cos(tt) * y * y + 2. * x * x * y) / 4.);
			}
			else
				return double(0.);

		// Rectangular fdtd::sources::signal_value(Top Hat) Analytical Derivative
		case envelope::rectangular:
			return double(0.);

		// Continuous Wave Analytical Derivative
		case envelope::gaussian_periods: // parameter1: the fundamental frequency
			period = double(1 / parameter1);
			t_0 = double(0.5) * parameter2 * period;
			delta_t = double(0.1) * parameter2 * period;

			if ((t > double(0.)) && (t < parameter2 * period))
			{
				exponent = exp(-(t - t_0) * (t - t_0) / (delta_t * delta_t));
				return double(- 2. * ePeak * (t - t_0) * exponent / (delta_t * delta_t));
			}
			else
				return double(0.);

		// Long Sine Wave Analytical Derivative
		case envelope::arctan_ramp:
			period = double(1 / parameter1);
			t_0 = 0.5 * parameter2 * period;
			return double(ePeak * period / (2. * ((t - t_0) * (t - t_0) + period * period)));

		// Gaussian Envelope Analytical Derivative
		case envelope::gaussian:
			exponent = exp(-2.7726 * (double(t - parameter1) / double(parameter2)) * (double(t - parameter1) / double(parameter2)));
			return double(- 2. * 2.7726 * ePeak * double(t - parameter1) * exponent / double(parameter2 * parameter2));

		// Undefined
		case envelope::linear_ramp:
			if (t < parameter2)
				return ePeak / parameter2 * parameter1;
			else
				return 0.;

		// Default -- No envelope
		default:
			return double(1.);
	}
}

/**
 * @brief the derivatives of all the Carrier functions
 * 
 * @param type [in] integer value corresponding to carrier type
 * @param t [in] time the carrier is being evaluated at
 * @param wavelength [in] wavelength of the carrier
 * @param offset [in] time/phase delay for the carrier
 * @return double 
 */
double d_carriers(int type, const double& t, const double& wavelength, const double& offset) {
	double carrier;
	switch (type)
	{
		case carrier::cosine:
			carrier = - 2. * PI * wavelength * sin(2. * PI * wavelength * (t - offset)); //cos(2.0 * PI * wavelength * t - offset);
			return carrier;
		
		case carrier::sine:
			carrier = wavelength * cos(wavelength * (t - offset));//sin(wavelength * t - offset);
			return carrier;

		default:
			return (double(0.));
	}
}

/**
 * @brief the derivative of the Complex Carrier function
 */ 
std::complex<double> d_complex_carriers(int type, const double& t, const double& wavelength, const double& chirp1, const double& chirp2, const std::complex<double>& iUnit, const double& delay, const double& width) {
	std::complex<double> carrier;
	std::complex<double> exponent1;
	std::complex<double> exponent2;
	std::complex<double> exponent3;
	switch (type)
	{
	case carrier::chirped:
		exponent1 = std::complex<double>(exp(-iUnit * 2.0 * PI * wavelength * t));
		exponent2 = std::complex<double>(exp(iUnit * chirp1 * 2.7726 * (double(t - delay) / double(width)) * (double(t - delay) / double(width))));
		exponent3 = std::complex<double>(exp(iUnit * chirp2 * 2.7726 * (double(t - delay) / double(width)) * (double(t - delay) / double(width)) * (double(t - delay) / double(width))));
		
		carrier = std::complex<double>(iUnit * exponent1 * exponent2 * exponent3 * (- 2 * PI * wavelength +
		2 * chirp1 * 2.7726 * double(t - delay) / double(width * width) +
		3 * chirp2 * 2.7726 * double(t - delay) * double(t - delay) / double(width * width * width)));
		return carrier;
	default:
		return std::complex<double>(1.);
	}
}

/**
 * @brief Analytic derivatives of signal function, Electric Field Source
 * 
 * @param t time in seconds
 * @return double derivative signal value in V/m
 * 
 * @details This function uses global vars to determine signal at a given time
 * 't'. This function handles: \n
 * - Raised Cosine \n
 * - Derived Raised Cosine \n
 * - CW Pulse \n
 * - Continuous Wave \n
 * - Modulated Gaussian \n
 * - Constant field \n
 * 
 */
double fdtd::sources::Signal::derivative(double t) const
{
	double period;
	double k_period;

	if (_shape == Shape::raised_cosine)
	{
		return (d_envelopes(envelope::raised_cosine,t,_k, _period,_e_peak));
	}
	//derived raised cosine
	else if (_shape == Shape::raised_cosine_derivative)
	{
		return (d_envelopes(envelope::raised_cosine_derivative, t, _k, _period, _e_peak));
		
	}
	else if (_shape == Shape::rectangular)
	{
		return (d_envelopes(envelope::rectangular, t, _dt, _period, _e_peak));
		
	}
	//cw pulse for SHG or THG
	else if (_shape == Shape::cw_pulse)
	{
		period = double(1 / _fundamental); // was left uninitialised here, roadmap §7
		k_period = double(PI2 / period);
		return (d_envelopes(envelope::gaussian_periods, t, _fundamental, _cw_periods, _e_peak) * carriers(carrier::sine, t, k_period, double(double(0.5) * _cw_periods * period))
		+ envelopes(envelope::gaussian_periods, t, _fundamental, _cw_periods, _e_peak) * d_carriers(carrier::sine, t, k_period, double(double(0.5) * _cw_periods * period)));
	}
	//long sin wave
	else if (_shape == Shape::continuous_wave)
	{
		k_period = double(PI2 * _fundamental);
		return (d_envelopes(envelope::arctan_ramp, t, _fundamental, _start_cw, _e_peak) * carriers(carrier::sine, t, k_period, double(0))
		+ envelopes(envelope::arctan_ramp, t, _fundamental, _start_cw, _e_peak) * d_carriers(carrier::sine, t, k_period, double(0)));
		
	}
	else if (_shape == Shape::gaussian)
	{
			if(_chirp == 0 && _chirp_2 == 0)
			{
				double pulse;
				pulse = (d_envelopes(envelope::gaussian, t, _delay, _width, _e_peak) * carriers(carrier::cosine,t,_f_0,double(0)) 
				+ envelopes(envelope::gaussian, t, _delay, _width, _e_peak) * d_carriers(carrier::cosine,t,_f_0,double(0)));
				return pulse;
			}
			else
			{
				std::complex<double> pulse;
				pulse = std::complex<double>( d_envelopes(envelope::gaussian, t, _delay, _width, _e_peak) * complex_carriers(carrier::chirped, t, _f_0,_chirp,_chirp_2,I_UNIT,_delay,_width ) 
				+ envelopes(envelope::gaussian, t, _delay, _width, _e_peak) * d_complex_carriers(carrier::chirped, t, _f_0,_chirp,_chirp_2,I_UNIT,_delay,_width )) ;
				return pulse.real();
			}

	}
	else if (_shape == Shape::gaussian_alt)
	{
		double pulse;
		double pulse2;
		if (t < 3 * (_delay + 2 * _width))
		{
			pulse = (d_envelopes(envelope::gaussian, t, _delay, _width, _e_peak) * carriers(carrier::cosine, t, _f_0, double(0))
			+ envelopes(envelope::gaussian, t, _delay, _width, _e_peak) * d_carriers(carrier::cosine, t, _f_0, double(0)));
			pulse2 = (d_envelopes(envelope::gaussian, t, _delay, _width, _e_peak) * carriers(carrier::cosine, t, _f_02, double(0))
			+ envelopes(envelope::gaussian, t, _delay, _width, _e_peak) * d_carriers(carrier::cosine, t, _f_02, double(0)));
		}
		else
		{
			pulse = 0.0;
			pulse2 = 0.0;
		}

		return pulse + pulse2;
	}

	else if (_shape == Shape::constant)
	{
		double pulse;
		if (_shape == Shape::constant)
			pulse = d_envelopes(envelope::linear_ramp,t,_dt, 1000.0 * _dt,_e_peak);

		return pulse;
	}
	else
		return double(0.);
}

/// Samples the signal at steps 0..steps: its peak, and the last step above threshold * peak.
fdtd::sources::SignalExtent fdtd::sources::Signal::extent(int steps, double threshold) const
{
	const double dt = _dt;
	SignalExtent e{0, 0, 0};
	for (int n = 0; n <= steps; ++n)
	{
		const double s = std::abs(value(n * dt));
		if (s > e.peak)
			e.peak = s, e.peak_step = n;
	}
	for (int n = steps; n >= 0; --n)
		if (std::abs(value(n * dt)) > threshold * e.peak)
		{
			e.last_on = n;
			break;
		}
	return e;
}

/**
 * @brief Whether a run of `steps` time steps covers the signal
 *
 * It fits when it has decayed to below 1e-3 of its peak before the last 2 % of the run. Steady
 * signals, continuous wave, constant field, always fit.
 * @param why  set to the reason when the signal does not fit
 */
bool fdtd::sources::Signal::fits(int steps, std::string &why) const
{
	if (steady())
		return true;
	const SignalExtent e = extent(steps);
	if (e.peak == 0)
	{
		why = "the source signal is zero over the " + std::to_string(steps) + " time steps of the run";
		return false;
	}
	if (e.last_on > steps - std::max(1, steps / 50))
	{
		why = "the source signal is still on at the end of the run (" + std::to_string(steps) + " time steps, signal above 1e-3 of its peak at step " +
			  std::to_string(e.last_on) + "): increase [grid] time_steps (see fdtd signal), or set [source] allow_truncated_signal = true";
		return false;
	}
	return true;
}

fdtd::sources::Signal::Signal(const Params &p)
{
	Signal &g = *this;
	g._shape = p.shape;
	g._dt = p.dt;
	g._e_peak = p.e_peak;
	g._fundamental = p.fundamental;
	g._start_cw = p.start_cw;
	g._chirp = p.chirp;
	g._chirp_2 = p.chirp_2;
	g._cw_periods = p.cw_periods;
	switch (g._shape)
	{
	case Shape::raised_cosine: g._period = double(1.) / p.max_frequency; break;
	case Shape::raised_cosine_derivative: g._period = double(2.) / p.max_frequency; break;
	case Shape::rectangular: g._period = double(1.) / p.max_frequency; break;
	case Shape::gaussian:
	case Shape::gaussian_alt:
		g._width = p.pulse_width;
		g._delay = p.pulse_width > 40e-15 ? 2.0 * g._width : 3 * g._width;
		g._f_0 = C0 / (p.center_wavelength);
		break;
	case Shape::sinc:
		/* See Emadeldeen Hassan and Antonio Calà Lesina, "Topology optimization of dispersive plasmonic
		   nanostructures in the time-domain," Opt. Express 30, 19557-19572, 2022 */
		g._period = 25 * 1e-15; // broadband
		break;
	default: break;
	}
	g._k = PI2 / g._period;
	if (g._shape == Shape::continuous_wave)
		g._cw_periods = 100;
}
