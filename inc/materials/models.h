/**
 * @file materials/models.h
 * @brief The dispersive material models and the dispatch from a material index to its model.
 *
 * Adding a model: a header in materials/models/ with the shape of models/cell.h, a Kind and its keys in
 * materials/definition.h, its parameters in load_slots, and one line in with_model.
 */
#ifndef FDTD_MATERIALS_MODELS_H
#define FDTD_MATERIALS_MODELS_H

#include "materials/slots.h"
#include "materials/models/dielectric.h"
#include "materials/models/drude.h"
#include "materials/models/drude_2cp.h"
#include "materials/models/drude_lorentz.h"
#include "materials/models/hydrodynamic.h"
#include "materials/models/lorentz.h"
#include "materials/models/magnetic_lorentz.h"
#include "materials/models/modified_lorentz.h"

namespace fdtd::materials
{
/**
 * @brief Calls f(Model{}) with the model of dispersive slot mm; nothing for a dielectric
 *
 * Called once per stretch: f instantiates its loop for the model, so the cells see no dispatch.
 */
template <class F>
void with_model(int mm, F &&f)
{
	switch (kind(mm))
	{
	case Kind::drude_cp: f(models::Drude2CP{}); break;
	case Kind::drude: f(models::Drude{}); break;
	case Kind::lorentz: f(models::Lorentz{}); break;
	case Kind::modified_lorentz: f(models::ModifiedLorentz{}); break;
	case Kind::hydrodynamic: f(models::Hydrodynamic{}); break;
	case Kind::drude_lorentz: f(models::DrudeLorentz{}); break;
	case Kind::magnetic_lorentz: f(models::MagneticLorentz{}); break;
	case Kind::dielectric: break;
	}
}

/// a Drude or Lorentz slot: its E update leaves the CPML and TF/SF terms to be added after it, scaled by ce2
inline bool adds_terms_after(int mm) { return is(mm, Kind::drude) || is(mm, Kind::lorentz); }

/// ce2 of a Drude or Lorentz slot, adds_terms_after
inline double after_terms_ce2(int mm)
{
	return is(mm, Kind::drude) ? coeffs<models::Drude>(mm).ce2 : coeffs<models::Lorentz>(mm).ce2;
}

/// the curl factor of slot m's E update, dt / eps: ce_2 of a dielectric, ce2 of a Drude or Lorentz slot;
/// false for the models without a plain factor
inline bool e_factor(int m, double &c)
{
	if (!dispersive(m))
	{
		c = model_parameters().ce_2[m];
		return true;
	}
	if (!adds_terms_after(m))
		return false;
	c = after_terms_ce2(m);
	return true;
}

/// the electron-gas coefficients of slot mm, hydrodynamic or Drude-Lorentz; zero for the other models
inline const ElectronGas &electron_gas(int mm)
{
	static const ElectronGas none{};
	if (is(mm, Kind::hydrodynamic))
		return coeffs<models::Hydrodynamic>(mm);
	if (is(mm, Kind::drude_lorentz))
		return coeffs<models::DrudeLorentz>(mm);
	return none;
}
} // namespace fdtd::materials

#endif
