/**
 * @file materials/thermal/thermal.h
 * @brief The two-temperature model: electron and lattice temperatures driven by the absorbed power,
 *        and the temperature-dependent Drude parameters they feed back, [temperature] model
 *
 * Known issue: update_ttm crashes, roadmap §7; the model awaits its review.
 */
#ifndef FDTD_MATERIALS_THERMAL_THERMAL_H
#define FDTD_MATERIALS_THERMAL_THERMAL_H

#include "materials/slots.h"

namespace fdtd::materials::thermal
{
/// the material the two-temperature model heats: the first Drude-Lorentz slot, else slot 0
inline int heated_slot()
{
	const int s = first_slot(Kind::drude_lorentz);
	return s < 0 ? 0 : s;
}

/// setup.cc
void load_probe_temperature_1d();
void load_static_temperature();
/// two_temperature.cc: the electron and lattice temperatures and the temperature-dependent material constants
bool update_ttm(int step); ///< step: the time step, pump-probe bookkeeping
void update_ttm_rk();
double electron_phonon_coupling(double);
double electron_heat_capacity(double);
double lattice_heat_capacity(double);
double plasma_frequency(double);
double damping_constant(double, double);
/// output.cc
void read_temperature();
void save_temperature(int step);
} // namespace fdtd::materials::thermal

#endif
