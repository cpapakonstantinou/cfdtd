/**
 * @file monitors/exact_multipoles.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @date 2026-07-20
 *
 * @brief Exact, beyond long-wavelength, multipole decomposition of the induced
 * current density inside the antenna volume.
 *
 * @details Implements the exact multipole moments of Alaee, Rockstuhl &
 * Fernandez-Corbaton, Opt. Commun. 407, 17-21, 2018, Table 2:
 *
 *   ED:  p_a   = -1/(iw) { Int J_a j0(kr) + (k^2/2) Int [3(r.J)r_a - r^2 J_a] j2(kr)/(kr)^2 }
 *   MD:  m_a   = (3/2) Int (r x J)_a j1(kr)/(kr)
 *   EQ:  Qe_ab = -3/(iw) { Int [3(r_b J_a + r_a J_b) - 2(r.J)d_ab] j1(kr)/(kr)
 *                        + 2k^2 Int [5 r_a r_b (r.J) - (r_a J_b + r_b J_a) r^2 - r^2 (r.J) d_ab] j3(kr)/(kr)^3 }
 *   MQ:  Qm_ab = 15 Int { r_a (r x J)_b + r_b (r x J)_a } j2(kr)/(kr)^2
 *
 * The induced, contrast, current is built from the volume E-field DFT:
 *
 *   J_w(r) = -i w eps0 [eps_r(r,w) - eps_host(w)] E_w(r)
 *
 * using the code's e^{+iwt} DFT convention, Re += E cos(wt), Im += E sin(wt),
 * for which d/dt -> -iw. Substituting J = -iw P_eq with
 * P_eq = eps0 (eps_r - eps_host) E_w, the -1/(iw) prefactors cancel and the
 * moments are evaluated directly in terms of P_eq, no division by omega.
 *
 * k in the spherical Bessel arguments is the host wavenumber
 * k = Re[n(first_medium, w)] * w / c0. The multipole origin is the nominal
 * antenna center, center_position_x, center_position_y - antenna_shift/dy,
 * center_position_z, matching the geometry-creation convention.
 *
 * Requires the volume E-field DFT, VolumeDft::re/im[E], which is accumulated when
 * dfg_optim || save_output || ABS_SIGMA || exact_multipoles. Works with or
 * without the "Advanced Model" flag: there is NO dependence on the material
 * model, hydrodynamics, TTM, etc.
 *
 * SHG/THG (dfg_optim = 3..6): the upper half of analysis_omega[] holds the harmonic
 * band, and the same loop yields the harmonic multipoles from the harmonic
 * internal field, "E-contrast" moments. In addition, the instantaneous
 * nonlinear polarization P_NL(r,t) = eps0 chi2:EE, SHG, or eps0 chi3:EEE, THG
 * is Fourier-transformed on the fly (sample_nonlinear_polarisation) and a second set
 * of moments is produced from the TOTAL current
 * J_tot = -iw [P_eq + P_NL,w], written to a separate file, so both options
 * are available for comparison.
 *
 * Long-wavelength, LWA, moments, [multipoles] lwa: the same integrals in the limit kr -> 0
 * j0 -> 1, j1/x -> 1/3, j2/x^2 -> 1/15, j3/x^3 -> 1/105 and the k^2 terms vanish, i.e. Alaee et al.
 * Table 1: p = Int P, m = (1/2) Int r x J, Qe_ab = Int [3(r_a P_b + r_b P_a) - 2 (r.P) d_ab],
 * Qm_ab = Int [r_a (r x J)_b + r_b (r x J)_a]. Like the exact set they are built from the volume
 * E-field DFT and the material permittivities, so they work with every material model.
 *
 * Output, root process, reports_dir, same style as cross_sections.txt:
 *   multipoles_exact.txt                 all frequencies, E-contrast current
 *   multipoles_exact_full_current.txt    harmonic band only, contrast + nl polarization
 *   multipoles_lwa.txt, multipoles_lwa_full_current.txt   the LWA moments
 *   <each>_cross_sections.txt            per-multipole scattering cross sections
 *
 * Column layout, 37 columns:
 *   LAMBDA[nm]  Re/Im px py pz  Re/Im mx my mz
 *   Re/Im Qe_xx yy zz xy xz yz  Re/Im Qm_xx yy zz xy xz yz
 *
 * Units: p [C m], m [A m^2], Qe [C m^2], Qm [A m^3], all per the DFT
 * normalization in use, fundamental band: normalized to the incident-pulse
 * spectrum, i.e. "per unit incident field"; harmonic band: raw DFT, identical
 * to the convention of the ART / Q-coefficient outputs.
 *
 * @copyright Copyright 2026
 */

#include "monitors/exact_multipoles.h"
#include "io/text.h"
#include "monitors/state.h"
#include "monitors/volume_dft.h"
#include "materials/material_map.h"
#include "materials/nonlinear/state.h"
#include "grid/state.h"
#include "sources/state.h"
#include "materials/materials.h"
#include "materials/nonlinear/nonlinear.h"
#include "grid/subdomain.h"
#include <string>
#include <vector>
#include "core/limits.h"
#include "grid/lattice.h"
#include "monitors/spec.h"
using fdtd::constants::C0;
using fdtd::constants::EP0;
using fdtd::constants::I_UNIT;
using fdtd::constants::PI2;
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;
using fdtd::materials::matind_x;
using fdtd::materials::matind_y;
using fdtd::materials::matind_z;

namespace
{
const fdtd::monitors::Spec &sp = fdtd::monitors::spec();        ///< what the monitors need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::monitors; // transitional: module state, monitors/state.h

namespace fdtd::monitors
{

// Spherical bessel weight functions

/**
 * @brief Computes j0(x), j1(x)/x, j2(x)/x^2, j3(x)/x^3 in one shot.
 *
 * @details The ratios, rather than the bare Bessel functions, are what enter
 * the Table-2 integrands, and they are finite at x->0, where they reduce to
 * the long-wavelength limits 1, 1/3, 1/15, 1/105, recovering Table 1.
 * A Taylor series is used for small arguments to avoid catastrophic
 * cancellation in the closed forms.
 */
static inline void sph_bessel_weights(const double &x, double &j0v, double &f1, double &f2, double &f3)
{
	if (x < double(0.05))
	{
		double x2 = x * x;
		double x4 = x2 * x2;
		j0v = double(1.) - x2 / double(6.) + x4 / double(120.);
		f1 = double(1.) / double(3.) - x2 / double(30.) + x4 / double(840.);
		f2 = double(1.) / double(15.) - x2 / double(210.) + x4 / double(7560.);
		f3 = double(1.) / double(105.) - x2 / double(1890.) + x4 / double(83160.);
	}
	else
	{
		double sx = double(sin(x));
		double cx = double(cos(x));
		double inv = double(1.) / x;
		double inv2 = inv * inv;
		double inv3 = inv2 * inv;
		double inv4 = inv2 * inv2;

		j0v = sx * inv;
		double j1 = sx * inv2 - cx * inv;
		double j2 = (double(3.) * inv3 - inv) * sx - double(3.) * inv2 * cx;
		double j3 = (double(15.) * inv4 - double(6.) * inv2) * sx - (double(15.) * inv3 - inv) * cx;

		f1 = j1 * inv;
		f2 = j2 * inv2;
		f3 = j3 * inv3;
	}
	return;
}

// Material helpers

/**
 * @brief Decides whether a cell/component belongs to the multipole
 * integration volume.
 *
 * @details Default, linear runs, dfg_optim < 3: the antenna material only.
 * For SHG/THG runs (dfg_optim = 3..6) the gap material and the shell are
 * treated as part of the emitting volume, since the bulk around a plasmonic
 * structure can act as the source of the nonlinear emission.
 */
static inline bool include_in_multipole_volume(const int &mat)
{
	const int base = mat;

	if (base == sp.antenna_material)
		return true;

	if (sp.shg || sp.thg) // dfg_optim == 3,4, THG, or 5,6, SHG
	{
		if (base == sp.gap_material)
			return true;
		if (base == sp.shell_material)
			return true;
	}
	return false;
}

/**
 * @brief The slot of a nonlinear material, whose chi tensor the nonlinear updates use,
 *        nonlinear::update_polarisation and the exact dielectric update; -1 for a linear one
 */
static inline int nl_chi_index(const int &mat)
{
	return fdtd::materials::nonlinear::nl[mat] ? mat : -1;
}

// Running dft of the nonlinear polarization

/**
 * @brief Additive DFT step for the instantaneous nonlinear polarization
 * P_NL(r,t), harmonic band only (m >= nr_an_frequencies/2).
 *
 * @details Called once per step by the P_NL monitor, monitors/monitor.cc, when
 * exact_multipoles && (SHG || THG). The polarization is evaluated straight
 * from the current E fields with the same chi2/chi3 contraction and the same
 * sta_ir_casing / PER-COMPONENT branching as the nonlinear update itself
 * nonlinear::polarisation_term, inc/materials/nonlinear/nonlinear.h, so the DFT'd
 * quantity is exactly the P_NL the simulation injects.
 *
 * @param t current simulation time of the electric field
 */
void sample_nonlinear_polarisation(double t)
{
	VolumeDft &V = volume_dft();
	if (!(sp.shg || sp.thg) || sp.nr_frequencies <= 2)
		return;

	int i;
	int j;
	int k;
	int I;
	int J;
	int K;
	int m;
	int mx;
	int my;
	int mz;
	int cx;
	int cy;
	int cz;
	double pnlx;
	double pnly;
	double pnlz;
	double s;

	const int mstart = sp.nr_frequencies / 2;

	// The harmonic-band P_NL DFT arrays must have been allocated for m>=mstart
	// setupSimulation.cpp. If they are still NULL, the per-timestep write
	// below would segfault at the first step; instead emit a clear one-time
	// message and skip the P_NL DFT, the contrast-current multipoles are
	// unaffected. This means the allocation block of RePNL*/ImPNL* was not
	// added, or its gate does not match this run.
	if (V.re_pnl[0] == nullptr || V.re_pnl[1] == nullptr || V.re_pnl[2] == nullptr ||
		V.im_pnl[0] == nullptr || V.im_pnl[1] == nullptr || V.im_pnl[2] == nullptr)
	{
		static bool warned = false;
		if (!warned && fdtd::grid::world_rank == ROOTPROCESS)
		{
			fdtd::io::log("monitors", "[Exact_Multipoles] WARNING: RePNL*/ImPNL* not allocated for the "
				   "harmonic band (m>=%d of %d). Skipping the nonlinear-polarization DFT; "
				   "the full-current output file will not be produced. Add the RePNL*/ImPNL* "
				   "allocation in setupSimulation.cpp.\n",
				   mstart, sp.nr_frequencies);
			fflush(stdout);
			warned = true;
		}
		return;
	}

	// precompute the DFT kernels once per timestep
	static double ucos[fdtd::max_frequencies];
	static double vsin[fdtd::max_frequencies];
	for (m = mstart; m < sp.nr_frequencies; ++m)
	{
		s = double(fdtd::sources::analysis_omega[m] * t);
		ucos[m] = double(cos(s));
		vsin[m] = double(sin(s));
	}

	for (i = 0; i < subdomain.inner_dim_x; ++i)
	{
		I = i + subdomain.inner_org_x;
		for (j = 0; j < subdomain.inner_dim_y; ++j)
		{
			J = j + subdomain.inner_org_y;
			for (k = 0; k < subdomain.inner_dim_z; ++k)
			{
				K = k + subdomain.inner_org_z;

				mx = int(matind_x[i][j][k]);
				my = int(matind_y[i][j][k]);
				mz = int(matind_z[i][j][k]);

				cx = nl_chi_index(mx);
				cy = nl_chi_index(my);
				cz = nl_chi_index(mz);

				if (cx < 0 && cy < 0 && cz < 0)
					continue;

				pnlx = double(0.);
				pnly = double(0.);
				pnlz = double(0.);

				if (cx >= 0)
				{
					pnlx = materials::nonlinear::polarisation_term(cx, I, J, K, 0, mx, my, mz, false);
				}
				if (cy >= 0)
				{
					pnly = materials::nonlinear::polarisation_term(cy, I, J, K, 1, mx, my, mz, false);
				}
				if (cz >= 0)
				{
					pnlz = materials::nonlinear::polarisation_term(cz, I, J, K, 2, mx, my, mz, false);
				}

				for (m = mstart; m < sp.nr_frequencies; ++m)
				{
					V.re_pnl[0][m][i][j][k] += pnlx * ucos[m];
					V.re_pnl[1][m][i][j][k] += pnly * ucos[m];
					V.re_pnl[2][m][i][j][k] += pnlz * ucos[m];

					V.im_pnl[0][m][i][j][k] += pnlx * vsin[m];
					V.im_pnl[1][m][i][j][k] += pnly * vsin[m];
					V.im_pnl[2][m][i][j][k] += pnlz * vsin[m];
				}
			}
		}
	}
	return;
}

// Moment accumulation kernel

// number of complex moments per frequency: p(3) + m(3) + Qe(6) + Qm(6)
static const int NR_MOMENTS = 18;

/**
 * @brief Adds the contribution of one cell to the 18 complex moments.
 *
 * @param P    equivalent polarization vector P = -J/(iw) at the cell [C/m^2]
 * @param r    cell position relative to the multipole origin [m]
 * @param kmag host wavenumber at this frequency [1/m]
 * @param miw  -i*analysis_omega[m], converts P back to J for the magnetic moments
 * @param acc  36 REALs: interleaved Re/Im of
 *             [px py pz mx my mz Qe(xx yy zz xy xz yz) Qm(xx yy zz xy xz yz)]
 */
static inline void accumulate_moments(const std::complex<double> P[3], const double r[3], const double &kmag,
									 const std::complex<double> &miw, double *acc)
{
	static const int PAIR_A[6] = {0, 1, 2, 0, 0, 1};
	static const int PAIR_B[6] = {0, 1, 2, 1, 2, 2};

	double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
	double rr = double(sqrt(r2));
	double j0v;
	double f1;
	double f2;
	double f3;
	sph_bessel_weights(kmag * rr, j0v, f1, f2, f3);

	double k2 = kmag * kmag;

	std::complex<double> rP = r[0] * P[0] + r[1] * P[1] + r[2] * P[2];

	std::complex<double> cross[3];
	cross[0] = r[1] * P[2] - r[2] * P[1];
	cross[1] = r[2] * P[0] - r[0] * P[2];
	cross[2] = r[0] * P[1] - r[1] * P[0];

	int a;
	int b;
	int q;
	std::complex<double> val;

	// electric dipole (T2-1):  p_a = Int [ P_a j0 + (k^2/2)(3(r.P)r_a - r^2 P_a) j2/(kr)^2 ]
	for (a = 0; a < 3; ++a)
	{
		val = P[a] * j0v + (double(0.5) * k2) * (double(3.) * rP * r[a] - r2 * P[a]) * f2;
		acc[2 * a] += val.real();
		acc[2 * a + 1] += val.imag();
	}

	// magnetic dipole (T2-2):  m_a = (3/2) Int (r x J)_a j1/(kr),  J = miw * P
	for (a = 0; a < 3; ++a)
	{
		val = double(1.5) * miw * cross[a] * f1;
		acc[2 * (3 + a)] += val.real();
		acc[2 * (3 + a) + 1] += val.imag();
	}

	// electric quadrupole (T2-3):
	// Qe_ab = 3 Int [ (3(r_b P_a + r_a P_b) - 2(r.P)d_ab) j1/(kr)
	//               + 2k^2 (5 r_a r_b (r.P) - (r_a P_b + r_b P_a) r^2 - r^2 (r.P) d_ab) j3/(kr)^3 ]
	for (q = 0; q < 6; ++q)
	{
		a = PAIR_A[q];
		b = PAIR_B[q];
		double dab = (a == b) ? double(1.) : double(0.);
		std::complex<double> sym = r[b] * P[a] + r[a] * P[b];
		val = double(3.) * ((double(3.) * sym - double(2.) * rP * dab) * f1 +
						  double(2.) * k2 * (double(5.) * r[a] * r[b] * rP - sym * r2 - r2 * rP * dab) * f3);
		acc[2 * (6 + q)] += val.real();
		acc[2 * (6 + q) + 1] += val.imag();
	}

	// magnetic quadrupole (T2-4):
	// Qm_ab = 15 Int { r_a (r x J)_b + r_b (r x J)_a } j2/(kr)^2,  J = miw * P
	for (q = 0; q < 6; ++q)
	{
		a = PAIR_A[q];
		b = PAIR_B[q];
		val = double(15.) * miw * (r[a] * cross[b] + r[b] * cross[a]) * f2;
		acc[2 * (12 + q)] += val.real();
		acc[2 * (12 + q) + 1] += val.imag();
	}
	return;
}

// Multipole scattering cross sections

/**
 * @brief Per-multipole scattering cross sections from a 36-real moment row,
 * following Alaee et al., 2018, Eq. (1):
 *
 *   Cpref = k^4 / (6 pi eps0^2 |E0|^2)
 *   Cp    = Cpref |p|^2
 *   Cm    = Cpref |m|^2 / c^2
 *   Cpq   = Cpref (k^2/120) ( sum_diag |Qe|^2 + 2 sum_offdiag |Qe|^2 )
 *   Cmq   = Cpref (k^2/120) ( sum_diag |Qm|^2 + 2 sum_offdiag |Qm|^2 ) / c^2
 *
 * The stored moments are already normalized to the incident spectrum on the
 * fundamental band (dft_normalization), so |E0| = 1 in those units and Cp..Cmq
 * are the physical per-unit-incident-field scattering cross sections. On a
 * harmonic band, raw DFT, there is no incident field to normalize against, so
 * the same combination is a raw analog, not a literal cross section.
 *
 * Moment row layout, idx 0..17, each = Re,Im pair:
 *   p: 0,1,2   m: 3,4,5   Qe: 6..11 (xx yy zz xy xz yz)   Qm: 12..17
 */
static inline void multipole_cross_sections(const double *M, const double &k, double &Cp, double &Cm, double &Cpq, double &Cmq)
{
	const double E0amp = double(1.0); // moments are per-unit-incident-field, see note above
	double k2 = k * k;
	double k4 = k2 * k2;
	double Cpref = k4 / (double(3.) * PI2 * EP0 * EP0 * E0amp * E0amp); // 6 pi = 3 * (2 pi) = 3 * PI2
	double c2 = C0 * C0;

#define NRM2(idx) (M[2 * (idx)] * M[2 * (idx)] + M[2 * (idx) + 1] * M[2 * (idx) + 1])
	Cp = Cpref * (NRM2(0) + NRM2(1) + NRM2(2));
	Cm = Cpref * (NRM2(3) + NRM2(4) + NRM2(5)) / c2;
	Cpq = Cpref * (k2 / double(120.)) *
		  (NRM2(6) + NRM2(7) + NRM2(8) + double(2.) * (NRM2(9) + NRM2(10) + NRM2(11)));
	Cmq = Cpref * (k2 / double(120.)) *
		  (NRM2(12) + NRM2(13) + NRM2(14) + double(2.) * (NRM2(15) + NRM2(16) + NRM2(17))) / c2;
#undef NRM2
	return;
}

/**
 * @brief Writes a cross-section file: LAMBDA, Cp, Cm, Cpq, Cmq, Ctot and the
 * multipole fractions, one row per frequency in [m0, m1).
 */
static void write_cross_section_file(const std::string &path, const double *M, int m0, int m1, int STRIDE)
{
	FILE *f = fdtd::io::open_report(path, "LAMBDA[nm]\tCp\tCm\tCpq\tCmq\tCtot\tf_pd\tf_md\tf_pq\tf_mq");
	if (f == NULL)
		return;
	for (int m = m0; m < m1; ++m)
	{
		double k = double(fdtd::sources::analysis_omega[m] / C0 * fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real());
		double Cp;
		double Cm;
		double Cpq;
		double Cmq;
		multipole_cross_sections(M + m * STRIDE, k, Cp, Cm, Cpq, Cmq);
		double Ctot = Cp + Cm + Cpq + Cmq;
		double inv = (Ctot != double(0.)) ? double(1.) / Ctot : double(0.);
		fprintf(f, "% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\n",
				(PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), Cp, Cm, Cpq, Cmq, Ctot,
				Cp * inv, Cm * inv, Cpq * inv, Cmq * inv);
	}
	fclose(f);
	return;
}

/// one row per frequency in [m0, m1): LAMBDA[nm] and the 36 moment values, Re, Im pairs
static void write_moment_file(const std::string &path, const double *M, int m0, int m1, int STRIDE)
{
	FILE *f = fdtd::io::open_report(path, "LAMBDA[nm]\t"
										  "Re_px\tIm_px\tRe_py\tIm_py\tRe_pz\tIm_pz\t"
										  "Re_mx\tIm_mx\tRe_my\tIm_my\tRe_mz\tIm_mz\t"
										  "Re_Qexx\tIm_Qexx\tRe_Qeyy\tIm_Qeyy\tRe_Qezz\tIm_Qezz\t"
										  "Re_Qexy\tIm_Qexy\tRe_Qexz\tIm_Qexz\tRe_Qeyz\tIm_Qeyz\t"
										  "Re_Qmxx\tIm_Qmxx\tRe_Qmyy\tIm_Qmyy\tRe_Qmzz\tIm_Qmzz\t"
										  "Re_Qmxy\tIm_Qmxy\tRe_Qmxz\tIm_Qmxz\tRe_Qmyz\tIm_Qmyz");
	if (f == NULL)
		return;
	for (int m = m0; m < m1; ++m)
	{
		fprintf(f, "% .8e", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9));
		for (int n = 0; n < STRIDE; ++n)
			fprintf(f, "\t% .8e", M[m * STRIDE + n]);
		fprintf(f, "\n");
	}
	fclose(f);
}

// Main entry point

/**
 * @brief Computes and writes the exact multipole moments.
 *
 * @details Call once after the time stepping is finished and AFTER the DFT normalisation, monitors' normalize
 * with the Q coefficients and ART spectra. All ranks must enter, collective
 * MPI_Reduce; ranks that own only CPML cells contribute zeros.
 */
void write_exact_multipoles()
{
	VolumeDft &V = volume_dft();
	if (!multipoles_wanted())
		return;
	const bool exact = exact_multipoles;
	const bool lwa = lwa_multipoles;

	int i;
	int j;
	int k;
	int m;

	const int STRIDE = 2 * NR_MOMENTS; // 36 REALs per frequency
	const int NTOT = sp.nr_frequencies * STRIDE;
	const int mstart = sp.nr_frequencies / 2;

	// The volume E-field DFT is mandatory. If it is NULL the allocation gate in
	// setupSimulation.cpp was not extended with ||exact_multipoles, integration
	// step 4. Report and bail rather than dereferencing NULL.
	if (!subdomain.cpml_only && V.re[0][0] == nullptr)
	{
		if (fdtd::grid::world_rank == ROOTPROCESS)
		{
			fdtd::io::log("monitors", "[Exact_Multipoles] ERROR: volume E-field DFT (ReTFEx) is not allocated. "
				   "Add ||EXACT_MULTIPOLES to the allocation gate in setupSimulation.cpp "
				   "(integration step 4) and to the accumulation gate in updateDFT.cpp "
				   "(step 5). No multipole output written.\n");
			fflush(stdout);
		}
		return;
	}

	// Include the nonlinear-polarization, full current, moments only if those
	// arrays were actually allocated; otherwise fall back to contrast-only.
	bool nl_run = (sp.shg || sp.thg) && sp.nr_frequencies > 2;
	if (nl_run && !subdomain.cpml_only && V.re_pnl[0] == nullptr)
	{
		if (fdtd::grid::world_rank == ROOTPROCESS)
		{
			fdtd::io::log("monitors", "[Exact_Multipoles] WARNING: RePNL*/ImPNL* not allocated; writing "
				   "contrast-current multipoles only (no full-current file).\n");
			fflush(stdout);
		}
		nl_run = false;
	}

	// one set of moment rows per kind, exact, LWA, and current, contrast, contrast + nl polarisation
	std::vector<double> loc(NTOT);
	std::vector<double> tot(NTOT);
	std::vector<double> locNL(nl_run ? NTOT : 0);
	std::vector<double> totNL(nl_run ? NTOT : 0);
	std::vector<double> locL(NTOT);
	std::vector<double> totL(NTOT);
	std::vector<double> locLNL(nl_run ? NTOT : 0);
	std::vector<double> totLNL(nl_run ? NTOT : 0);

	if (!subdomain.cpml_only)
	{
		const double dV = lat.d[0] * lat.d[1] * lat.d[2];

		// multipole origin = nominal antenna center used by the geometry builders
		const int cX = sp.centre[0];
		const int cY = sp.centre[1] - int(std::floor((sp.antenna_shift / lat.d[1]) + 0.5));
		const int cZ = sp.centre[2];

		double r[3];
		std::complex<double> P[3];
		std::complex<double> Pt[3];
		int mx;
		int my;
		int mz;
		bool inx;
		bool iny;
		bool inz;

		for (m = 0; m < sp.nr_frequencies; ++m)
		{
			const double kmag = double(fdtd::sources::analysis_omega[m] / C0 * fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real());
			const std::complex<double> eps_host = fdtd::materials::complex_permittivity(sp.medium_first, fdtd::sources::analysis_omega[m]);
			const std::complex<double> miw = -I_UNIT * fdtd::sources::analysis_omega[m];

			double *acc = loc.data() + m * STRIDE;
			double *accL = locL.data() + m * STRIDE;
			double *accNL = nl_run ? (locNL.data() + m * STRIDE) : NULL;
			double *accLNL = nl_run ? (locLNL.data() + m * STRIDE) : NULL;
			const bool do_nl = nl_run && (m >= mstart);

			for (i = 0; i < subdomain.inner_dim_x; ++i)
			{
				r[0] = double(i + subdomain.output_start[0] - cX) * lat.d[0];

				for (j = 0; j < subdomain.inner_dim_y; ++j)
				{
					r[1] = double(j + subdomain.output_start[1] - cY) * lat.d[1];

					for (k = 0; k < subdomain.inner_dim_z; ++k)
					{
						mx = int(matind_x[i][j][k]);
						my = int(matind_y[i][j][k]);
						mz = int(matind_z[i][j][k]);

						inx = include_in_multipole_volume(mx);
						iny = include_in_multipole_volume(my);
						inz = include_in_multipole_volume(mz);

						bool has_nl = do_nl && (V.re_pnl[0][m][i][j][k] != double(0.) || V.im_pnl[0][m][i][j][k] != double(0.) ||
												V.re_pnl[1][m][i][j][k] != double(0.) || V.im_pnl[1][m][i][j][k] != double(0.) ||
												V.re_pnl[2][m][i][j][k] != double(0.) || V.im_pnl[2][m][i][j][k] != double(0.));

						if (!inx && !iny && !inz && !has_nl)
							continue;

						r[2] = double(k + subdomain.output_start[2] - cZ) * lat.d[2];

						// equivalent polarization of the contrast current:
						// P_eq = eps0 (eps_r - eps_host) E_w, per Yee component
						P[0] = std::complex<double>(0.);
						P[1] = std::complex<double>(0.);
						P[2] = std::complex<double>(0.);

						if (inx)
							P[0] = EP0 * (fdtd::materials::complex_permittivity(mx, fdtd::sources::analysis_omega[m]) - eps_host) *
								   std::complex<double>(V.re[0][0][m][i][j][k], V.im[0][0][m][i][j][k]);
						if (iny)
							P[1] = EP0 * (fdtd::materials::complex_permittivity(my, fdtd::sources::analysis_omega[m]) - eps_host) *
								   std::complex<double>(V.re[0][1][m][i][j][k], V.im[0][1][m][i][j][k]);
						if (inz)
							P[2] = EP0 * (fdtd::materials::complex_permittivity(mz, fdtd::sources::analysis_omega[m]) - eps_host) *
								   std::complex<double>(V.re[0][2][m][i][j][k], V.im[0][2][m][i][j][k]);

						if (inx || iny || inz)
						{
							if (exact)
								accumulate_moments(P, r, kmag, miw, acc);
							if (lwa)
								accumulate_moments(P, r, double(0.), miw, accL); // the kr -> 0 limit
						}

						// total current on the harmonic band: contrast + nl polarization
						if (do_nl)
						{
							Pt[0] = P[0] + std::complex<double>(V.re_pnl[0][m][i][j][k], V.im_pnl[0][m][i][j][k]);
							Pt[1] = P[1] + std::complex<double>(V.re_pnl[1][m][i][j][k], V.im_pnl[1][m][i][j][k]);
							Pt[2] = P[2] + std::complex<double>(V.re_pnl[2][m][i][j][k], V.im_pnl[2][m][i][j][k]);
							if (exact)
								accumulate_moments(Pt, r, kmag, miw, accNL);
							if (lwa)
								accumulate_moments(Pt, r, double(0.), miw, accLNL);
						}
					}
				}
			}
		}

		const double dVall = dV;
		for (auto *v : {&loc, &locNL, &locL, &locLNL})
			for (double &x : *v)
				x *= dVall;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (exact)
	{
		MPI_Reduce(loc.data(), tot.data(), NTOT, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if (nl_run)
		{
			MPI_Reduce(locNL.data(), totNL.data(), NTOT, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	if (lwa)
	{
		MPI_Reduce(locL.data(), totL.data(), NTOT, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
		if (nl_run)
			MPI_Reduce(locLNL.data(), totLNL.data(), NTOT, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (fdtd::grid::world_rank == ROOTPROCESS)
	{
		const std::string dir = sp.reports_dir + "/";
		const int nf = sp.nr_frequencies;
		if (exact)
		{
			write_moment_file(dir + "multipoles_exact.txt", tot.data(), 0, nf, STRIDE);
			write_cross_section_file(dir + "multipoles_exact_cross_sections.txt", tot.data(), 0, nf, STRIDE); // Alaee et al. Eq. 1
			if (nl_run)
			{
				write_moment_file(dir + "multipoles_exact_full_current.txt", totNL.data(), mstart, nf, STRIDE);
				write_cross_section_file(dir + "multipoles_exact_full_current_cross_sections.txt", totNL.data(), mstart, nf, STRIDE);
			}
		}
		if (lwa)
		{
			write_moment_file(dir + "multipoles_lwa.txt", totL.data(), 0, nf, STRIDE);
			write_cross_section_file(dir + "multipoles_lwa_cross_sections.txt", totL.data(), 0, nf, STRIDE);
			if (nl_run)
			{
				write_moment_file(dir + "multipoles_lwa_full_current.txt", totLNL.data(), mstart, nf, STRIDE);
				write_cross_section_file(dir + "multipoles_lwa_full_current_cross_sections.txt", totLNL.data(), mstart, nf, STRIDE);
			}
		}
	}
	return;
}
} // namespace fdtd::monitors
