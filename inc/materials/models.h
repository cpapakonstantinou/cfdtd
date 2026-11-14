/**
 * @file materials/models.h
 * @brief The dispersive material models and the dispatch from a material index to its model.
 *
 * Adding a model: a header in materials/models/ with the shape of models/cell.h, a code range in
 * materials/codes.h, and one line in with_model.
 */
#ifndef FDTD_MATERIALS_MODELS_H
#define FDTD_MATERIALS_MODELS_H

#include "materials/codes.h"
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
 * @brief Calls f(Model{}) with the model of dispersive material index mm, code - DispMateOffsetD;
 *        nothing for an index past the last model
 *
 * Called once per stretch: f instantiates its loop for the model, so the cells see no dispatch.
 */
template <class F>
void with_model(int mm, F &&f)
{
	if (mm < DrudeCP)
		f(models::Drude2CP{});
	else if (mm < Drude)
		f(models::Drude{});
	else if (mm < Lorentz)
		f(models::Lorentz{});
	else if (mm < ModLorentz)
		f(models::ModifiedLorentz{});
	else if (mm < NonLinFreeElectron)
		f(models::Hydrodynamic{});
	else if (mm < NonlinearDrudeLorentz)
		f(models::DrudeLorentz{});
	else if (mm < MagneticLorentz)
		f(models::MagneticLorentz{});
}

/// the electron-gas coefficients of dispersive material index mm, hydrodynamic or nonlinear Drude-Lorentz;
/// zero for the other models
inline const ElectronGas &electron_gas(int mm)
{
	static const ElectronGas none{};
	if (mm >= ModLorentz && mm < NonLinFreeElectron)
		return coeffs<models::Hydrodynamic>(mm);
	if (mm >= NonLinFreeElectron && mm < NonlinearDrudeLorentz)
		return coeffs<models::DrudeLorentz>(mm);
	return none;
}
} // namespace fdtd::materials

#endif
