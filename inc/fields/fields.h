/** \file   fields/fields.h
 *  \brief  The electromagnetic fields of this rank and their time history.
 *
 *  Yee grid, staggered by half a cell: E components on cell edges, H components on cell faces.
 *  Every array is indexed [i][j][k] over the rank's cells plus the halo, 0 .. cells along each
 *  axis, k contiguous.
 *  @verbatim
      set linear    the fields of the run
      set harmonic  the nonlinear, harmonic, fields, when they are kept apart from the linear ones
                    [nonlinear] split; the hydrodynamic model
      E_prev        E at the previous step, per set: second-order ADEs, Drude + critical points,
      E_prev_2      E two steps back: magnetic Lorentz
      H_prev        H at the previous step: magnetic force term of the hydrodynamic model
      E_scratch     E - coeff dP_NL before it is committed, nonlinear E approach, or the copy of E
                    for the magneto-optic cross term, magnetic Lorentz; double, not double
   @endverbatim
 *  Until the update module owns its arguments, the legacy names (Ex .. Hz, ex_nl, ex_prev, ...) are
 *  references into this object.
 */
#ifndef FDTD_FIELDS_FIELDS_H
#define FDTD_FIELDS_FIELDS_H

#include "core/carray.h"

namespace fdtd::fields
{
/// array indices used by Fields
namespace index
{
constexpr int x = 0; ///< field component
constexpr int y = 1; ///< field component
constexpr int z = 2; ///< field component
constexpr int linear = 0; ///< field set
constexpr int harmonic = 1; ///< field set
} // namespace index

/// which optional arrays a run needs
struct Needs
{
	bool harmonic = false;	///< separate harmonic fields, and their E_prev
	bool e_prev_2 = false;	///< E two steps back
	bool h_prev = false;	///< H at the previous step
	bool e_scratch = false; ///< E_scratch
};

/// The fields of this rank
class Fields
{
public:
	ctensor<double> E[2][3];		  ///< [set][component], V/m
	ctensor<double> H[2][3];		  ///< [set][component], A/m
	ctensor<double> E_prev[2][3];  ///< [set][component]
	ctensor<double> E_prev_2[3];	  ///< [component], linear set
	ctensor<double> H_prev[3];	  ///< [component], linear set
	ctensor<double> E_scratch[3];	  ///< [component]

	/**
	 * @brief Allocates, zeroed, the arrays of a rank with cells (nx, ny, nz)
	 * @return false when out of memory
	 */
	bool allocate(int nx, int ny, int nz, const Needs &needs);

	/// frees every array
	void release();
};

/// this rank's fields, until the solver object owns them
Fields &fields();

/// The components by name, on this rank, the transitional names the update code uses
inline ctensor<double> &Ex = fields().E[0][0];                                   ///< Electric field x component on this rank, Yee edge, V/m, indexed [i][j][k] with halo.
inline ctensor<double> &Ey = fields().E[0][1];                                   ///< Electric field y component on this rank, Yee edge, V/m, indexed [i][j][k] with halo.
inline ctensor<double> &Ez = fields().E[0][2];                                   ///< Electric field z component on this rank, Yee edge, V/m, indexed [i][j][k] with halo.
inline ctensor<double> &Hx = fields().H[0][0];                                   ///< Magnetic field x component on this rank, Yee face, A/m, indexed [i][j][k] with halo.
inline ctensor<double> &Hy = fields().H[0][1];                                   ///< Magnetic field y component on this rank, Yee face, A/m, indexed [i][j][k] with halo.
inline ctensor<double> &Hz = fields().H[0][2];                                   ///< Magnetic field z component on this rank, Yee face, A/m, indexed [i][j][k] with halo.
inline ctensor<double> &ex_prev = fields().E_prev[0][0];                              ///< Ex at the previous step, needed by second-order ADEs.
inline ctensor<double> &ey_prev = fields().E_prev[0][1];                              ///< Ey at the previous step, needed by second-order ADEs.
inline ctensor<double> &ez_prev = fields().E_prev[0][2];                              ///< Ez at the previous step, needed by second-order ADEs.
inline ctensor<double> &ex_temp = fields().E_scratch[0];                               ///< Magnetic Lorentz scratch copy of Ex for the magneto-optic cross term.
inline ctensor<double> &ey_temp = fields().E_scratch[1];                               ///< Magnetic Lorentz scratch copy of Ey for the magneto-optic cross term.
inline ctensor<double> &ez_temp = fields().E_scratch[2];                               ///< Magnetic Lorentz scratch copy of Ez for the magneto-optic cross term.
inline ctensor<double> &ex_prev_2 = fields().E_prev_2[0];                            ///< Ex two steps back, magnetic Lorentz.
inline ctensor<double> &ey_prev_2 = fields().E_prev_2[1];                            ///< Ey two steps back, magnetic Lorentz.
inline ctensor<double> &ez_prev_2 = fields().E_prev_2[2];                            ///< Ez two steps back, magnetic Lorentz.
inline ctensor<double> &ex_nl = fields().E[1][0];                                ///< Harmonic, nonlinear, Ex field when linear and nonlinear fields are split.
inline ctensor<double> &ey_nl = fields().E[1][1];                                ///< Harmonic, nonlinear, Ey field when linear and nonlinear fields are split.
inline ctensor<double> &ez_nl = fields().E[1][2];                                ///< Harmonic, nonlinear, Ez field when linear and nonlinear fields are split.
inline ctensor<double> &hx_nl = fields().H[1][0];                                ///< Harmonic, nonlinear, Hx field when linear and nonlinear fields are split.
inline ctensor<double> &hy_nl = fields().H[1][1];                                ///< Harmonic, nonlinear, Hy field when linear and nonlinear fields are split.
inline ctensor<double> &hz_nl = fields().H[1][2];                                ///< Harmonic, nonlinear, Hz field when linear and nonlinear fields are split.
inline ctensor<double> &hx_prev = fields().H_prev[0];                              ///< Hx at the previous step, magnetic Lorentz force term.
inline ctensor<double> &hy_prev = fields().H_prev[1];                              ///< Hy at the previous step, magnetic Lorentz force term.
inline ctensor<double> &hz_prev = fields().H_prev[2];                              ///< Hz at the previous step, magnetic Lorentz force term.
inline ctensor<double> &ex_prev_nl = fields().E_prev[1][0];                           ///< Harmonic Ex at the previous step.
inline ctensor<double> &ey_prev_nl = fields().E_prev[1][1];                           ///< Harmonic Ey at the previous step.
inline ctensor<double> &ez_prev_nl = fields().E_prev[1][2];                           ///< Harmonic Ez at the previous step.
} // namespace fdtd::fields

#endif // FDTD_FIELDS_FIELDS_H
