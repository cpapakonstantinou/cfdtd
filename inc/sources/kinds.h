/**
 * @file sources/kinds.h
 * @brief The source codes: `[source] type`, `[source] signal`, and the envelopes and carriers the
 *        signals are built from, names as in the INI file where there is a key
 */
#ifndef FDTD_SOURCES_KINDS_H
#define FDTD_SOURCES_KINDS_H

namespace fdtd::sources
{
/// [source] type
namespace kind
{
inline constexpr int none = 0;
inline constexpr int plane_wave = 1;
inline constexpr int dipole = 2;
inline constexpr int gaussian_beam = 3;
inline constexpr int wave_superposition = 4;
inline constexpr int mode = 5; ///< waveguide mode on a plane, sources/mode_source.cc
} // namespace kind

/// [source] signal
namespace pulse
{
inline constexpr int raised_cosine = 1;
inline constexpr int raised_cosine_derivative = 2;
inline constexpr int continuous_wave = 3;
inline constexpr int cw_pulse = 4;
inline constexpr int rectangular = 5;
inline constexpr int gaussian = 6;
inline constexpr int gaussian_alt = 7;
inline constexpr int constant = 8;
inline constexpr int sinc = 11;
} // namespace pulse

/// envelope of a signal
namespace envelope
{
inline constexpr int raised_cosine = 1;            ///< (1 - cos)^3 / 8 over one period
inline constexpr int raised_cosine_derivative = 2; ///< its derivative
inline constexpr int rectangular = 3;              ///< top hat
inline constexpr int gaussian_periods = 4;         ///< Gaussian over a number of carrier periods, CW pulse
inline constexpr int arctan_ramp = 5;              ///< slow turn-on to a constant amplitude, continuous wave
inline constexpr int gaussian = 6;                 ///< Gaussian of a given width, FWHM, and delay
inline constexpr int linear_ramp = 8;              ///< linear build-up to a constant field
} // namespace envelope

/// carrier of a signal
namespace carrier
{
inline constexpr int cosine = 1;  ///< cos(2 pi f, t - offset)
inline constexpr int sine = 2;    ///< sin(k, t - offset)
inline constexpr int chirped = 1; ///< the complex carrier with linear and quadratic chirp
} // namespace carrier
} // namespace fdtd::sources

#endif
