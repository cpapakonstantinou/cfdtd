/** \file   src/sources/field_scale.cc
 *  \brief  The scale of each frequency of the field outputs, from the source spectrum.
 */
#include "sources/field_scale.h"
#include "core/constants.h"
#include "sources/spec.h"
#include "sources/state.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();
} // namespace

std::string fdtd::sources::FieldScale::description() const
{
	char text[240];
	if (order == 1)
		std::snprintf(text, sizeof text, "E over the incident field at %.4g nm%s", wavelength * 1e9, phase ? "" : ", amplitudes only");
	else
		std::snprintf(text, sizeof text, "E at %.4g nm, harmonic order %d, over the peak incident amplitude, at %.4g nm%s", wavelength * 1e9, order,
					  reference_wavelength * 1e9, phase ? ", pulse delay taken off the phase" : ", amplitudes only");
	return text;
}

h5io::FieldNote fdtd::sources::FieldScale::note() const
{
	h5io::FieldNote n;
	n.scale = factor;
	n.wavelength = wavelength;
	n.normalization = normalization();
	n.order = order;
	n.reference_wavelength = reference_wavelength;
	n.reference_amplitude = reference_amplitude;
	n.reference_phase = reference_phase;
	n.dt = sp.dt;
	n.description = description();
	return n;
}

fdtd::sources::FieldScale fdtd::sources::field_scale(int m, int set, bool phase, bool stored_normalised)
{
	using fdtd::constants::C0;
	using fdtd::constants::PI2;
	// the pump band is the first fundamental_point_width.size() points, the harmonic band follows
	const int pumps = int(fundamental_point_width.size());
	FieldScale s;
	s.phase = phase;
	s.pump = m < pumps || pumps == 0 ? m : m - pumps;
	const int band_order = s.pump == m ? 1 : int(std::lround(analysis_omega[m] / analysis_omega[s.pump]));
	s.order = band_order * (set == 1 ? 2 : 1);
	s.wavelength = PI2 * C0 / ((set == 1 ? 2 : 1) * analysis_omega[m]);

	const std::complex<double> pump = dft_normalization[s.pump];
	std::complex<double> reference;
	if (s.order == 1)
	{
		reference = phase ? pump : std::complex<double>(std::abs(pump));
		s.reference_wavelength = PI2 * C0 / analysis_omega[s.pump];
	}
	else
	{
		// the largest pump amplitude, the phase of this point's pump to the power k
		int peak = 0;
		for (int j = 1; j < std::max(pumps, 1); ++j)
			if (amplitude_normalization[j] > amplitude_normalization[peak])
				peak = j;
		reference = std::polar(amplitude_normalization[peak], phase ? s.order * std::arg(pump) : 0.0);
		s.reference_wavelength = PI2 * C0 / analysis_omega[peak];
	}
	s.reference_amplitude = std::abs(reference) * sp.dt;
	s.reference_phase = std::arg(reference);

	// the pump band, stored normalised, keeps its values exactly
	if (s.order == 1 && stored_normalised)
		return s;
	const std::complex<double> stored = !stored_normalised ? 1.0 : phase ? dft_normalization[m] : std::complex<double>(amplitude_normalization[m]);
	s.factor = stored / reference;
	return s;
}
