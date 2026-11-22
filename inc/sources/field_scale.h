/**
 * @file sources/field_scale.h
 * @brief How the field outputs present each frequency: relative to the incident field
 *
 * Every frequency of an output is k times a pump point j of the [spectrum] band: k = 1 for the pump
 * band, the harmonic order for the harmonic band, twice that for the hydrodynamic harmonic set.
 *
 *   k = 1:  E / E_inc(w_j)                                  the incident field at the same frequency
 *   k > 1:  E / (|E_inc|_peak exp(i k arg E_inc(w_j)))     a fraction of the peak incident amplitude
 *
 * E_inc is the DFT of the source signal, so both are in V s/m and the ratio has no unit. A harmonic is
 * not the response to one pump frequency, so the harmonic band takes one amplitude, the largest of the
 * pump points, and keeps its spectral shape. The phase k arg E_inc(w_j) is the one a harmonic inherits
 * from where the pulse sits on the run's clock, the delay of the pulse: dividing it off is what the pump
 * band gets from E_inc(w_j). The same holds with or without a nonlinearity, so a linear run with a
 * harmonic band shows the pump's own content at k w on the scale of a harmonic. With [output]
 * phase_normalize off, only amplitudes are divided.
 *
 * The spectral amplitude in V s/m is the value times reference_amplitude exp(i reference_phase); the
 * raw DFT, the sum over the steps of E exp(i w t), is that over dt.
 */
#ifndef FDTD_SOURCES_FIELD_SCALE_H
#define FDTD_SOURCES_FIELD_SCALE_H

#include "io/h5io.h"

#include <complex>
#include <string>

namespace fdtd::sources
{
struct FieldScale
{
	std::complex<double> factor = 1; ///< multiplies the stored DFT into the output
	int order = 1;                   ///< k: the output frequency over the pump point's
	int pump = 0;                    ///< j: the pump point
	double wavelength = 0;           ///< of the output [m]
	double reference_wavelength = 0; ///< of the pump point whose amplitude divides [m]
	double reference_amplitude = 0;  ///< the |E_inc| that divides [V s/m]
	double reference_phase = 0;      ///< the phase taken off [rad]
	bool phase = true;               ///< the phase is relative to the incident field

	/// "incident_field" or "pump_field"
	const char *normalization() const { return order == 1 ? "incident_field" : "pump_field"; }
	/// one line for the reader of a file
	std::string description() const;
	/// the scale and attributes of an HDF5 field file
	h5io::FieldNote note() const;
};

/// The scale of frequency point m of field set `set`, 0 the fields, 1 the hydrodynamic harmonic set at
/// twice the frequency. `stored_normalised`: the stored DFT is already divided by the normalisation of
/// point m, dft_normalization or amplitude_normalization, as the planes and the volume are.
FieldScale field_scale(int m, int set, bool phase, bool stored_normalised);
} // namespace fdtd::sources

#endif // FDTD_SOURCES_FIELD_SCALE_H
