/**
 * @file app/signal_config.cc
 * @brief The run's source signal from the configuration, app/signal_config.h
 */
#include "app/signal_config.h"
#include "app/run_config.h"
#include "monitors/state.h"
#include "sources/kinds.h"
#include "sources/signal.h"

void fdtd::app::configure_signal()
{
	namespace pulse = fdtd::sources::pulse;
	using Shape = fdtd::sources::Signal::Shape;
	RunConfig &rc = edit_run_config(); // configure phase
	const int type = rc.source.signal_type;
	rc.source.continuous_wave = type == pulse::continuous_wave;
	rc.source.cw_pulse = type == pulse::cw_pulse;

	fdtd::sources::Signal::Params p;
	p.shape = type == pulse::raised_cosine ? Shape::raised_cosine
			: type == pulse::raised_cosine_derivative ? Shape::raised_cosine_derivative
			: type == pulse::continuous_wave ? Shape::continuous_wave
			: type == pulse::cw_pulse ? Shape::cw_pulse
			: type == pulse::rectangular ? Shape::rectangular
			: type == pulse::gaussian ? Shape::gaussian
			: type == pulse::gaussian_alt ? Shape::gaussian_alt
			: type == pulse::constant ? Shape::constant
			: type == pulse::sinc ? Shape::sinc
								  : Shape::none;
	p.dt = rc.time.dt;
	p.e_peak = rc.source.e_peak;
	p.fundamental = rc.source.fundamental_frequency;
	p.start_cw = rc.source.start_cw;
	p.chirp = rc.source.chirp;
	p.chirp_2 = rc.source.chirp_2;
	p.cw_periods = rc.source.cw_periods;
	p.max_frequency = fdtd::monitors::max_source_frequency;
	p.pulse_width = rc.source.pulse_width;
	p.center_wavelength = rc.source.center_freq;
	fdtd::sources::set_run_signal(fdtd::sources::Signal(p));
}
