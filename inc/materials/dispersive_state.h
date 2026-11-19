/**
 * @file materials/dispersive_state.h
 * @brief The material models' state and parameters as two objects
 *
 * DispersiveState: what the ADE models carry from step to step, per E component: the polarisation P
 * and its history, the critical-point, pole, polarisations, the currents J, Drude-type models, the
 * hydrodynamic and two-temperature extensions, for the linear and the separate harmonic fields.
 * ModelParameters: the parameters of each material slot, load_slots fills them: permittivities, conductivities,
 * Drude and critical-point constants. One owner each, so allocation, release, restart and the
 * memory report walk one list (for_each_field, for_each_pole) instead of naming every array.
 *
 * The kernels and the set-up read the members (`mp.epre_inf[mm]`, `ds.c[0].p`) through a local
 * reference to model_parameters() / dispersive_state(). Only the hydrodynamic model's own arrays
 * still have per-component names, hydro::px_nl, hydro::grad_div_x, ...; materials/hydro/state.h,
 * until its x/y/z copies are folded.
 */
#ifndef FDTD_MATERIALS_DISPERSIVE_STATE_H
#define FDTD_MATERIALS_DISPERSIVE_STATE_H

#include "core/array.h"
#include "materials/slots.h"
#include "materials/material_field.h"

namespace fdtd::materials
{
/// the state of one E component, by material slot
struct ComponentState
{
	using Field = material_field<double, max_slots>;
	Field p; ///< P and its history, p_prev_2/3: magnetic Lorentz
	Field p_prev; ///< P and its history, p_prev_2/3: magnetic Lorentz
	Field p_prev_2; ///< P and its history, p_prev_2/3: magnetic Lorentz
	Field p_prev_3; ///< P and its history, p_prev_2/3: magnetic Lorentz
	Field j; ///< current density, linear, harmonic
	Field j_nl; ///< current density, linear, harmonic
	Field j_prev; ///< hydrodynamic current history and extra term
	Field j_extra; ///< hydrodynamic current history and extra term
	array3<double> p_nl_prev;           ///< harmonic P one step back, one array for all materials
	array4<double> pc[max_slots]; ///< critical-point poles and one step back
	array4<double> pc_prev[max_slots]; ///< critical-point poles and one step back
	array4<double> pc_nl[max_slots]; ///< the same for the harmonic fields
	array4<double> pc_prev_nl[max_slots]; ///< the same for the harmonic fields
	// the nonlocal and hydrodynamic models, all materials in one padded array, [i+1][j+1][k+1]
	array3<double> pc_sum; ///< sum over the poles, or P for hydrodynamics
	array3<double> pc_sum_nl; ///< sum over the poles, or P for hydrodynamics
	array3<double> p_nonlocal; ///< nonlocal P and its second-pass copy
	array3<double> p_nonlocal_2; ///< nonlocal P and its second-pass copy
	array3<double> grad_div; ///< grad(div P) now and two steps back, GNOR; unpadded
	array3<double> grad_div_n; ///< grad(div P) now and two steps back, GNOR; unpadded
	array3<double> grad_div_n_1; ///< grad(div P) now and two steps back, GNOR; unpadded
	array3<double> e_padded; ///< E copied into the padded layout
	array3<double> e_padded_nl; ///< E copied into the padded layout
	array3<double> electron_density;             ///< hydrodynamic electron density n along this component
	array3<double> eps_inf;                      ///< background permittivity per cell, hydrodynamic, padded

	/// every per-material field, in allocation order
	template <class F>
	void for_each_field(F f)
	{
		for (Field *x : {&j, &j_nl, &j_extra, &j_prev, &p, &p_prev, &p_prev_2, &p_prev_3})
			f(*x);
	}
	/// the linear pole arrays of material slot mm
	template <class F>
	void for_each_pole(int mm, F f)
	{
		f(pc[mm]), f(pc_prev[mm]);
	}
	/// frees every array
	void release()
	{
		for_each_field([](Field &x) { x.reset(); });
		for (int mm = 0; mm < max_slots; ++mm)
			for (array4<double> *x : {&pc[mm], &pc_prev[mm], &pc_nl[mm], &pc_prev_nl[mm]})
				x->reset();
		for (array3<double> *x : {&p_nl_prev, &pc_sum, &pc_sum_nl, &p_nonlocal, &p_nonlocal_2, &grad_div, &grad_div_n, &grad_div_n_1,
								   &e_padded, &e_padded_nl, &electron_density, &eps_inf})
			x->reset();
	}
};

struct DispersiveState
{
	ComponentState c[3]; ///< Ex, Ey, Ez
	template <class F>
	void for_each_field(F f)
	{
		for (ComponentState &s : c)
			s.for_each_field(f);
	}
	void release()
	{
		for (ComponentState &s : c)
			s.release();
	}
	/// what a restart file holds, in file order, each kind for x, y, z: P, its history, the harmonic
	/// P one step back, the poles, linear and harmonic, now and one step back, J, linear, harmonic
	template <class F>
	void for_each_saved(F f)
	{
#define FDTD_EACH(m) f(c[0].m), f(c[1].m), f(c[2].m)
		FDTD_EACH(p), FDTD_EACH(p_prev), FDTD_EACH(p_prev_2), FDTD_EACH(p_prev_3), FDTD_EACH(p_nl_prev);
		FDTD_EACH(pc), FDTD_EACH(pc_prev), FDTD_EACH(pc_nl), FDTD_EACH(pc_prev_nl);
		FDTD_EACH(j), FDTD_EACH(j_nl);
#undef FDTD_EACH
	}
	/// the linear pole arrays of slot mm, component by component, x, y, z; pole, then one step back
	template <class F>
	void for_each_pole(int mm, F f)
	{
		for (ComponentState &s : c)
			f(s.pc[mm]);
		for (ComponentState &s : c)
			f(s.pc_prev[mm]);
	}
};

/// the parameters of every material slot, materials/slots.h
struct ModelParameters
{
	double epre[max_slots];         ///< relative permittivity, dielectrics
	double sig_m[max_slots];        ///< conductivity [S/m]
	double ce_1[max_slots];         ///< E update coefficient of E^n, lossy dielectric
	double ce_2[max_slots];         ///< E update coefficient of curl H
	double n_refr_index[max_slots]; ///< real refractive index per slot, source delays
	double nr_cp[max_slots];         ///< critical points per slot
	double epre_inf[max_slots]; ///< eps_inf, static eps
	double epre_sta[max_slots]; ///< eps_inf, static eps
	double ome_plasma[max_slots]; ///< Drude plasma frequency and damping [rad/s]
	double gam_plasma[max_slots]; ///< Drude plasma frequency and damping [rad/s]
	double a_3[max_slots];           ///< magnetic Lorentz coupling term
	double big_omega_cp[max_slots][max_poles]; ///< pole frequency and damping
	double big_gamma_cp[max_slots][max_poles]; ///< pole frequency and damping
	double a_coeff_cp[max_slots][max_poles]; ///< pole amplitude and phase
	double phi_cp[max_slots][max_poles]; ///< pole amplitude and phase
	double plasma_freq_add_term;         ///< ITO: plasma frequency minus its two-temperature value at rest
	// the nonlocal and hydrodynamic models
	double non_loc[max_slots]; ///< nonlocal coefficient, beta, GNOR diffusion
	double beta_nl[max_slots]; ///< nonlocal coefficient, beta, GNOR diffusion
	double diff_coeff[max_slots]; ///< nonlocal coefficient, beta, GNOR diffusion
	double gnor_1[max_slots]; ///< GNOR update coefficients
	double gnor_2[max_slots]; ///< GNOR update coefficients
	double equilibrium_electron_density[max_slots];
	double effective_electron_mass[max_slots];
	double powerterm[max_slots];
	double inverse_edens[max_slots];
	double inverse_chargedens[max_slots];
};

extern DispersiveState the_dispersive_state; ///< this rank's model state
extern ModelParameters the_model_parameters; ///< the run's model parameters
inline DispersiveState &dispersive_state() { return the_dispersive_state; }
inline ModelParameters &model_parameters() { return the_model_parameters; }
} // namespace fdtd::materials

#endif
