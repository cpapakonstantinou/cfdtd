/**
 * @file materials/dispersive_state.h
 * @brief The material models' state and parameters as two objects
 *
 * DispersiveState: what the ADE models carry from step to step, per E component: the polarisation P
 * and its history, the critical-point, pole, polarisations, the currents J, Drude-type models, the
 * hydrodynamic and two-temperature extensions, for the linear and the separate harmonic fields.
 * ModelParameters: the per-code parameters set_parameters fills, permittivities, conductivities,
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

#include "core/carray.h"
#include "materials/codes.h"
#include "materials/material_field.h"

namespace fdtd::materials
{
/// the state of one E component, by material slot
struct ComponentState
{
	using Field = material_field<double, nr_disp_mat_d>;
	Field p; ///< P and its history, p_prev_2/3: magnetic Lorentz
	Field p_prev; ///< P and its history, p_prev_2/3: magnetic Lorentz
	Field p_prev_2; ///< P and its history, p_prev_2/3: magnetic Lorentz
	Field p_prev_3; ///< P and its history, p_prev_2/3: magnetic Lorentz
	Field j; ///< current density, linear, harmonic
	Field j_nl; ///< current density, linear, harmonic
	Field j_prev; ///< hydrodynamic current history and extra term
	Field j_extra; ///< hydrodynamic current history and extra term
	ctensor<double> p_nl_prev;           ///< harmonic P one step back, one array for all materials
	ctetrad<double> pc[nr_disp_mat_d]; ///< critical-point poles and one step back
	ctetrad<double> pc_prev[nr_disp_mat_d]; ///< critical-point poles and one step back
	ctetrad<double> pc_nl[nr_disp_mat_d]; ///< the same for the harmonic fields
	ctetrad<double> pc_prev_nl[nr_disp_mat_d]; ///< the same for the harmonic fields
	// the nonlocal and hydrodynamic models, all materials in one padded array, [i+1][j+1][k+1]
	ctensor<double> pc_sum; ///< sum over the poles, or P for hydrodynamics
	ctensor<double> pc_sum_nl; ///< sum over the poles, or P for hydrodynamics
	ctensor<double> p_nonlocal; ///< nonlocal P and its second-pass copy
	ctensor<double> p_nonlocal_2; ///< nonlocal P and its second-pass copy
	ctensor<double> grad_div; ///< grad(div P) now and two steps back, GNOR; unpadded
	ctensor<double> grad_div_n; ///< grad(div P) now and two steps back, GNOR; unpadded
	ctensor<double> grad_div_n_1; ///< grad(div P) now and two steps back, GNOR; unpadded
	ctensor<double> e_padded; ///< E copied into the padded layout
	ctensor<double> e_padded_nl; ///< E copied into the padded layout
	ctensor<double> electron_density;             ///< hydrodynamic electron density n along this component
	ctensor<double> eps_inf;                      ///< background permittivity per cell, hydrodynamic, padded

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
		for (int mm = 0; mm < nr_disp_mat_d; ++mm)
			for (ctetrad<double> *x : {&pc[mm], &pc_prev[mm], &pc_nl[mm], &pc_prev_nl[mm]})
				x->reset();
		for (ctensor<double> *x : {&p_nl_prev, &pc_sum, &pc_sum_nl, &p_nonlocal, &p_nonlocal_2, &grad_div, &grad_div_n, &grad_div_n_1,
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

/// per material code, dielectric codes 0..NRMAXMATERIALS-1, and per dispersive slot (0..nr_disp_mat_d-1)
struct ModelParameters
{
	double epre[NRMAXMATERIALS];         ///< relative permittivity, dielectric codes
	double sig_m[NRMAXMATERIALS];        ///< conductivity [S/m]
	double ce_1[NRMAXMATERIALS];         ///< E update coefficient of E^n, lossy dielectric
	double ce_2[NRMAXMATERIALS];         ///< E update coefficient of curl H
	double n_refr[100];                  ///< refractive index per dielectric code, epre = n_refr^2
	double n_refr_index[NRMAXMATERIALS]; ///< real refractive index per code, source delays
	double nr_cp[nr_disp_mat_d];         ///< critical points per slot
	double epre_inf[nr_disp_mat_d]; ///< eps_inf, static eps
	double epre_sta[nr_disp_mat_d]; ///< eps_inf, static eps
	double ome_plasma[nr_disp_mat_d]; ///< Drude plasma frequency and damping [rad/s]
	double gam_plasma[nr_disp_mat_d]; ///< Drude plasma frequency and damping [rad/s]
	double a_3[nr_disp_mat_d];           ///< modified-Lorentz term
	double big_omega_cp[nr_disp_mat_d][NrMaxCP]; ///< pole frequency and damping
	double big_gamma_cp[nr_disp_mat_d][NrMaxCP]; ///< pole frequency and damping
	double a_coeff_cp[nr_disp_mat_d][NrMaxCP]; ///< pole amplitude and phase
	double phi_cp[nr_disp_mat_d][NrMaxCP]; ///< pole amplitude and phase
	double plasma_freq_add_term;         ///< ITO: plasma frequency minus its two-temperature value at rest
	// the nonlocal and hydrodynamic models, per slot / code
	double non_loc[nr_disp_mat_d]; ///< nonlocal coefficient, beta, GNOR diffusion
	double beta_nl[nr_disp_mat_d]; ///< nonlocal coefficient, beta, GNOR diffusion
	double diff_coeff[nr_disp_mat_d]; ///< nonlocal coefficient, beta, GNOR diffusion
	double gnor_1[nr_disp_mat_d]; ///< GNOR update coefficients
	double gnor_2[nr_disp_mat_d]; ///< GNOR update coefficients
	double equilibrium_electron_density[nr_disp_mat_d];
	double effective_electron_mass[nr_disp_mat_d];
	double powerterm[NRMAXMATERIALS];
	double inverse_edens[NRMAXMATERIALS];
	double inverse_chargedens[NRMAXMATERIALS];
};

extern DispersiveState the_dispersive_state; ///< this rank's model state
extern ModelParameters the_model_parameters; ///< the run's model parameters
inline DispersiveState &dispersive_state() { return the_dispersive_state; }
inline ModelParameters &model_parameters() { return the_model_parameters; }
} // namespace fdtd::materials

#endif
