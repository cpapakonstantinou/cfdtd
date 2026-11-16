/**
 * @file materials/component.h
 * @brief The arrays of one E component, for the material models.
 *
 * A model is written once, for component A (0 x, 1 y, 2 z); it reaches its arrays through
 * Component<A> and the other two components through the cyclic axes B1 = A+1 and B2 = A+2.
 * Example: the magnetisation cross term of the magnetic Lorentz model is
 *   M_{B1} dE_{B2} - M_{B2} dE_{B1}, Ex: My dEz - Mz dEy, Ey: Mz dEx - Mx dEz, Ez: Mx dEy - My dEx
 *
 * The arrays are module state, inc/materials/state.h, inc/fields/fields.h; this header is where they are
 * named per component, so the models do not know the x/y/z names.
 */
#ifndef FDTD_MATERIALS_COMPONENT_H
#define FDTD_MATERIALS_COMPONENT_H

#include "grid/subdomain.h"
#include "core/views.h"
#include "core/carray.h"
#include "fields/fields.h"
#include "materials/hydro/state.h"
#include "materials/state.h"
#include "materials/dispersive_state.h"
#include "materials/spec.h"
#include "grid/lattice.h"

namespace fdtd::materials
{
/// x, y or z of three objects, chosen at compile time
template <int A, class T>
constexpr T &pick(T &x, T &y, T &z)
{
	static_assert(A >= 0 && A < 3, "axis is 0, 1 or 2");
	if constexpr (A == 0)
		return x;
	else if constexpr (A == 1)
		return y;
	else
		return z;
}

/// 1 / d along axis b
inline double inv_d(int b) { return b == 0 ? fdtd::grid::lattice().inv_d[0] : b == 1 ? fdtd::grid::lattice().inv_d[1] : fdtd::grid::lattice().inv_d[2]; }

/**
 * @brief The arrays of the E component A that the material models read and write
 *
 * Index conventions: e, p and the critical-point polarisations are indexed [i][j][k], with the
 * material first for p, pc; the arrays marked "padded" are indexed [i+1][j+1][k+1].
 */
template <int A>
struct Component
{
	static constexpr int B1 = (A + 1) % 3; ///< the other two axes, cyclic
	static constexpr int B2 = (A + 2) % 3; ///< the other two axes, cyclic
	static ComponentState &state() { return the_dispersive_state.c[A]; } ///< this component's model state

	// field and its history
	static auto &e() { return pick<A>(fdtd::fields::Ex, fdtd::fields::Ey, fdtd::fields::Ez); }
	static auto &e_prev() { return pick<A>(fdtd::fields::ex_prev, fdtd::fields::ey_prev, fdtd::fields::ez_prev); }
	static auto &e_prev_2() { return pick<A>(fdtd::fields::ex_prev_2, fdtd::fields::ey_prev_2, fdtd::fields::ez_prev_2); }
	static auto &e_scratch() { return pick<A>(fdtd::fields::ex_temp, fdtd::fields::ey_temp, fdtd::fields::ez_temp); } ///< copy before this step's update, magnetic Lorentz
	static auto &e_nl() { return pick<A>(fdtd::fields::ex_nl, fdtd::fields::ey_nl, fdtd::fields::ez_nl); }            ///< harmonic field
	static auto &e_prev_nl() { return pick<A>(fdtd::fields::ex_prev_nl, fdtd::fields::ey_prev_nl, fdtd::fields::ez_prev_nl); }

	// polarisation of the Drude / Lorentz pole, per material
	static auto &p() { return state().p; }
	static auto &p_prev() { return state().p_prev; }
	static auto &p_prev_2() { return state().p_prev_2; }
	static auto &p_prev_3() { return state().p_prev_3; }

	// critical-point, or Lorentz-pole, polarisations [material][pole]
	static auto &pc() { return state().pc; }
	static auto &pc_prev() { return state().pc_prev; }
	static auto &pc_nl() { return state().pc_nl; }
	static auto &pc_prev_nl() { return state().pc_prev_nl; }
	static auto &pc_sum() { return state().pc_sum; }          ///< sum over poles, or P for hydrodynamics, padded
	static auto &pc_sum_nl() { return state().pc_sum_nl; } ///< padded

	// nonlocal model
	static auto &p_nonlocal() { return state().p_nonlocal; } ///< padded
	static auto &p_nonlocal_prev() { return state().p_nl_prev; }
	static auto &p_nonlocal_2() { return state().p_nonlocal_2; } ///< padded
	static auto &grad_div() { return state().grad_div; }
	static auto &grad_div_n() { return state().grad_div_n; }
	static auto &grad_div_n_1() { return state().grad_div_n_1; }

	// free-electron current, hydrodynamic and nonlinear Drude-Lorentz models
	static auto &j() { return state().j; }
	static auto &j_nl() { return state().j_nl; }

	static auto &magnetization() { return pick<A>(fdtd::materials::magnetization_x, fdtd::materials::magnetization_y, fdtd::materials::magnetization_z); }
};

/// the value of v at (x[0], x[1], x[2]) shifted by da along axis a and db along axis b
inline double shifted(const view3<double> &v, const int (&x)[3], int a, int da, int b = 0, int db = 0)
{
	int y[3] = {x[0], x[1], x[2]};
	y[a] += da;
	y[b] += db;
	return v[y[0]][y[1]][y[2]];
}
} // namespace fdtd::materials

#endif
