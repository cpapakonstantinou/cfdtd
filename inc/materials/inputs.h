/**
 * @file materials/inputs.h
 * @brief What set_parameters needs from the run: the media, the user values of the built-in models,
 *        the custom material and the nonlinear tensors of the four media
 *
 * Built by the application from the configuration, app/setup.cc; the materials module reads it
 * instead of the run configuration, docs/dev/Roadmap.md §2.5.
 */
#ifndef FDTD_MATERIALS_INPUTS_H
#define FDTD_MATERIALS_INPUTS_H

namespace fdtd::materials
{
/// [custom_material]: the user's Drude + two critical-point material
struct CustomMaterial
{
	double eps_inf = 0;
	double plasma_frequency = 0;
	double plasma_damping = 0;
	double cp_frequency = 0;
	double cp_frequency_2 = 0;
	double cp_damping = 0;
	double cp_damping_2 = 0;
	double cp_amplitude = 0;
	double cp_amplitude_2 = 0;
	double cp_phase = 0;
	double cp_phase_2 = 0;
	double electron_density = 0;
	double fermi_velocity = 0;
};

struct Inputs
{
	int antenna = 0; ///< material codes of the structure and of the three layers
	int first = 0; ///< material codes of the structure and of the three layers
	int middle = 0; ///< material codes of the structure and of the three layers
	int second = 0; ///< material codes of the structure and of the three layers
	double n[4] = {};                 ///< refractive indices of the dielectric slots 0..3
	bool graphene_optim = false;      ///< materials 1, 2 are the averaging media of the graphene optimisation
	double middle_layer = 0;          ///< middle layer thickness [m]: selects the ITO film model (350)
	double graphene_stretch = 0;      ///< graphene (102): factor on the plasma frequency and the CP amplitudes
	double fermi_velocity = 0;        ///< [m/s] hydrodynamic beta
	double gnor_diffusion = 0;        ///< [m^2/s] GNOR diffusion
	double carrier_density = 0; ///< [m^-3] ITO
	double carrier_density_back = 0; ///< [m^-3] ITO
	CustomMaterial custom;
	/// chi2, Voigt, 18, and chi3 (21) tables of the structure, the first, middle and second medium
	double chi2[4][18] = {};
	double chi3[4][21] = {};
};
} // namespace fdtd::materials

#endif
