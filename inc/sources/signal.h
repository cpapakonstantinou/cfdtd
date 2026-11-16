/**
 * @file sources/signal.h
 * @brief The time signal of the sources and its spectrum
 *
 * @verbatim
      Signal                  [source] signal with its timing; value(t), derivative(t), signal.cc
      run_signal()            the run's signal, which the sources hold
      transform_pulse()       the DFT of the signal at the analysis frequencies: the normalisation
                              of every DFT monitor (dft_normalization, amplitude_normalization)
      build_plane_waves       the superposition list from the [wave.*] sections (plane_waves.cc)
   @endverbatim
 */
#ifndef FDTD_SOURCES_SIGNAL_H
#define FDTD_SOURCES_SIGNAL_H

#include <string>

#include <vector>

struct PlaneWave;
struct WaveFrame;

namespace fdtd::sources
{
/// where the signal is on within `steps` time steps
struct SignalExtent
{
	double peak;   ///< largest |signal|
	int peak_step; ///< step of the peak
	int last_on;   ///< last step with |signal| above threshold * peak
};

/** The run's source signal: `[source] signal` with its timing, as one value the sources hold.
 *  Built from Params, the application fills them from the configuration, app::configure_signal;
 *  evaluated by value(t) and derivative(t). */
class Signal
{
public:
	enum class Shape { none, raised_cosine, raised_cosine_derivative, rectangular, cw_pulse, continuous_wave, gaussian, gaussian_alt, constant, sinc };

	/// what a signal is made from, [source] signal and its [signal.*] section, the band
	struct Params
	{
		Shape shape = Shape::none;
		double dt = 0; ///< time step [s], peak field [V/m]
		double e_peak = 0; ///< time step [s], peak field [V/m]
		double fundamental = 0; ///< CW signals: frequency [Hz] and settling periods
		double start_cw = 0; ///< CW signals: frequency [Hz] and settling periods
		double chirp = 0; ///< Gaussian: linear and quadratic chirp
		double chirp_2 = 0; ///< Gaussian: linear and quadratic chirp
		int cw_periods = 30;                  ///< CW pulse length [periods]
		double max_frequency = 0;             ///< highest analysed frequency [Hz]: the band-limited signals' length
		double pulse_width = 0;               ///< Gaussian width [s]
		double center_wavelength = 0;         ///< Gaussian carrier [m]
	};
	Signal() = default;
	explicit Signal(const Params &p); ///< sets the timing, period, width, delay, carrier, from the shape

	double value(double t) const;      ///< envelope x carrier at time t
	double derivative(double t) const; ///< its time derivative, the derivative-pulse option

	Shape shape() const { return _shape; }
	bool has_carrier() const { return _shape == Shape::gaussian || _shape == Shape::gaussian_alt; } ///< its own carrier (center_wavelength)
	bool steady() const { return _shape == Shape::continuous_wave || _shape == Shape::constant; }     ///< no end
	double period() const { return _period; }               ///< length of the band-limited and sinc signals [s]
	double carrier_frequency() const { return _f_0; }       ///< carrier of the Gaussian signals [Hz]
	int cw_periods() const { return _cw_periods; }          ///< length of the CW signals [periods]

	SignalExtent extent(int steps, double threshold = 1e-3) const; ///< sampled at steps 0..steps
	/// a run of `steps` time steps covers the signal, it has decayed at the end; `why` says why not
	bool fits(int steps, std::string &why) const;

private:
	Shape _shape = Shape::none;
	double _dt = 0; ///< time step, peak field
	double _e_peak = 0; ///< time step, peak field
	double _period = 0; ///< band-limited signals: length and 2 pi / length
	double _k = 0; ///< band-limited signals: length and 2 pi / length
	double _delay = 0; ///< Gaussian envelope: centre and width [s]
	double _width = 0; ///< Gaussian envelope: centre and width [s]
	double _f_0 = 0; ///< Gaussian carriers [Hz], f_02: the second pulse of gaussian_alt, 0
	double _f_02 = 0; ///< Gaussian carriers [Hz], f_02: the second pulse of gaussian_alt, 0
	double _chirp = 0; ///< linear and quadratic chirp
	double _chirp_2 = 0; ///< linear and quadratic chirp
	double _fundamental = 0; ///< CW signals: frequency and settling periods
	double _start_cw = 0; ///< CW signals: frequency and settling periods
	int _cw_periods = 30;
};

/// the run's signal, set by the application, app::configure_signal
const Signal &run_signal();
void set_run_signal(const Signal &s);

/// the run's signal at time t: run_signal().value(t), monitor normalisation, tools
inline double signal_value(double t) { return run_signal().value(t); }
inline double signal_derivative(double t) { return run_signal().derivative(t); }
/// the DFT of the signal at the analysis frequencies; returns the run's step count, `steps`, or the
/// CW length for a continuous wave
int transform_pulse(int steps);
inline bool signal_fits(int steps, std::string &why) { return run_signal().fits(steps, why); }
inline SignalExtent signal_extent(int steps, double threshold = 1e-3) { return run_signal().extent(steps, threshold); }
bool build_plane_waves(const std::vector<PlaneWave> &input, std::vector<PlaneWave> &waves, WaveFrame &frame);
} // namespace fdtd::sources

#endif
