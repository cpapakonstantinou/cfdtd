/** \file   boundary/cpml.h
 *  \brief  Convolutional PML, CFS-CPML: grading profiles, auxiliary psi fields and their updates.
 *
 *  The CPML stretches each coordinate u in its absorbing slabs as
 *  @verbatim
      s_u = kappa_u + sigma_u / (alpha_u + j omega eps0)
      d/du  ->  (1 / kappa_u) d/du + psi_u,     psi_u^{n+1} = b_u psi_u^n + c_u d/du
      b_u = exp(-(sigma_u / kappa_u + alpha_u) dt / eps0)
      c_u = sigma_u (b_u - 1) / (kappa_u (sigma_u + kappa_u alpha_u))      (alpha_max != 0)
   @endverbatim
 *  graded polynomially from the interface (m = N) to the outer end (m = 0):
 *  @verbatim
      sigma(m) = sigma_max ((N - m) / N)^exponent,   sigma_max = 0.8 (exponent + 1) / (377 du)
      kappa(m) = 1 + (kappa_max - 1) ((N - m) / N)^exponent
      alpha(m) = alpha_max (m / N)^alpha_exponent
   @endverbatim
 *  Integer positions in m are the E samples along u, half positions (m + 1/2) the H samples;
 *  the half table's last entry (m = N) repeats the integer one.
 *
 *  The module owns the profiles of every axis and face and the psi arrays of every field,
 *  component, derivative axis and face, for the linear fields and, when they are kept apart,
 *  the harmonic ones. The users read it through cpml(), cached once per file in the hot paths.
 *
 *  Measured reflection: test/cpml/check.py, dipole, oblique, and test/cpml/plane_wave.py.
 */
#ifndef FDTD_BOUNDARY_CPML_H
#define FDTD_BOUNDARY_CPML_H

#include "core/carray.h"
#include "grid/subdomain.h"
#include "grid/stretching.h"

namespace fdtd::boundary
{
/// array indices used by Cpml
namespace index
{
constexpr int x = 0; ///< axis / field component
constexpr int y = 1; ///< axis / field component
constexpr int z = 2; ///< axis / field component
constexpr int near = 0; ///< face of an axis: low / high global index
constexpr int far = 1; ///< face of an axis: low / high global index
constexpr int integer = 0; ///< grading position: E samples (m) / H samples (m + 1/2)
constexpr int half = 1; ///< grading position: E samples (m) / H samples (m + 1/2)
constexpr int e = 0; ///< field
constexpr int h = 1; ///< field
constexpr int linear = 0; ///< field set, harmonic: the separate nonlinear fields
constexpr int harmonic = 1; ///< field set, harmonic: the separate nonlinear fields
} // namespace index

/// grading of one axis and face, m = 0, outer end .. N, interface, by grading position
struct Profile
{
	// 1 / kappa is the coordinates', grid/stretching.h
	cvector<double> b[2];			///< psi decay b
	cvector<double> c[2];			///< curl weight c
};

/// grading constants, the same on every axis today
struct Grading
{
	double exponent = 3;		 ///< polynomial order of sigma and kappa
	double kappa_max = 5;		 ///< kappa at the outer end
	double alpha_max = 0.05;	 ///< alpha at the interface (S/m)
	double alpha_exponent = 1; ///< polynomial order of alpha
};

/// The CPML of this rank
class Cpml
{
public:
	int thickness[3] = {};	///< cells per face, by axis
	Profile profile[3][2];	///< [axis][face]
	/// psi[set][field][component][derivative axis][face]; allocated only where the derivative axis
	/// has a CPML face on this rank and differs from the component
	ctensor<double> psi[2][2][3][3][2];
	bool harmonic = false; ///< the separate harmonic fields have psi too, set by allocate
	/// the material the E faces scale their psi with: the first / second medium on the y faces, the y
	/// layer's medium, the 1D grid's, on the x and z faces, or each cell's own, cellwise
	struct Background
	{
		int near_y = 0; ///< material codes of the y faces
		int far_y = 0; ///< material codes of the y faces
		bool cellwise = false;     ///< each cell's material, the maps extend into the CPML: mode source
	} background;

	/**
	 * @brief Computes the grading tables of every axis and face
	 * @param cells  CPML thickness per axis in cells
	 * @param du     cell size per axis in m
	 * @param dt     time step in s
	 */
	void set_profiles(const int cells[3], const double du[3], double dt, const Grading &g = {});

	/**
	 * @brief Allocates, zeroed, the psi arrays of the faces this rank touches
	 * @param near, far  this rank has the near / far CPML face of each axis
	 * @param harmonic   also the psi of the separate harmonic fields
	 * @return false when out of memory
	 *
	 * A face's arrays span the rank's cells, minus the inner cells along the face's axis.
	 */
	bool allocate(const grid::Subdomain &s, const bool near[3], const bool far[3], bool harmonic);

	/// frees the profiles and the psi arrays
	void release();
};

/// this rank's CPML, until the solver object owns it
Cpml &cpml();

/// psi updates of the E faces of Ex, y faces, z faces, Ey (x, z) and Ez (x, y), after the E scan
void update_psi_ex();
void update_psi_ey();
void update_psi_ez();

/**
 * @brief CPML psi and TF/SF terms that a dispersive model adds to its new E component A, scaled by its
 *        own curl coefficient c; updates the psi of the cell
 * @param p, q  the curl derivatives: curl = p - q, p = dH_{A+2}/d(A+1), q = dH_{A+1}/d(A+2)
 * @verbatim
      Ex: y faces, + p = dHz/dy, z faces (- q = dHy/dz), TF/SF y planes
      Ey: x faces, - q = dHz/dx, z faces (+ p = dHx/dz), TF/SF z and x planes
      Ez: y faces, - q = dHx/dy, x faces (+ p = dHy/dx), TF/SF y planes
   @endverbatim
 */
template <int A>
void cpml_tfsf_terms(int i, int j, int k, double c, double p, double q, double &e_new);
} // namespace fdtd::boundary

#endif // FDTD_BOUNDARY_CPML_H
