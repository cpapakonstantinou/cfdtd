/**
 * @file materials/spec.h
 * @brief The run switches the material models read: time step, which electron / temperature model
 *        is on, the harmonic fields
 *
 * Set by the application once the run is configured, after the time step, app/setup.cc; the model
 * kernels, inc/materials/models, and the coefficients read it instead of the run configuration.
 */
#ifndef FDTD_MATERIALS_SPEC_H
#define FDTD_MATERIALS_SPEC_H

namespace fdtd::materials
{
struct ModelSpec
{
	double dt = 0;                 ///< time step [s]
	double fermi_velocity = 0;     ///< [m/s] hydrodynamic pressure term
	bool shg = false; ///< harmonic generation on, the slots with a tensor are nonlinear
	bool thg = false; ///< harmonic generation on, the slots with a tensor are nonlinear
	bool harmonic = false;         ///< separate harmonic fields (*_nl)
	bool cellwise = false;         ///< every cell's own material in the CPML, mode source
	bool nonlocal = false;         ///< nonlocal electron model, hydrodynamic or GNOR
	bool hydrodynamics = false;    ///< hydrodynamic electron gas
	bool gnor = false;             ///< GNOR diffusion
	bool two_temperature = false;  ///< two-temperature model
	bool load_temperature = false; ///< temperature read from file
	bool pump = false; ///< the pump-probe sequence of the two-temperature model
	bool probe = false; ///< the pump-probe sequence of the two-temperature model
	bool pump_probe = false; ///< the pump-probe sequence of the two-temperature model
	int antenna_material = 0; ///< the structure's and the background's material slots
	int medium_first = 0; ///< the structure's and the background's material slots
	int medium_second = 0; ///< the structure's and the background's material slots
	bool staircase = false;        ///< one material per cell, [grid] staircase
	int tfsf_y_shrink = 0;         ///< TF/SF box shrink along y [cells], the hydrodynamic density window
};
extern ModelSpec the_model_spec; ///< set by set_model_spec
inline const ModelSpec &model_spec() { return the_model_spec; }
inline void set_model_spec(const ModelSpec &s) { the_model_spec = s; }
} // namespace fdtd::materials

#endif
