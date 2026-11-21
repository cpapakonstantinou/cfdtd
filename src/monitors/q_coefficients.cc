#include "monitors/q_coefficients.h"
#include "io/text.h"
#include "monitors/state.h"
#include "monitors/art.h"
#include "grid/subdomain.h"
#include "core/views.h"
#include "grid/state.h"
#include "geometry/state.h"
#include "sources/state.h"
#include "materials/materials.h"
#include "monitors/flux_boxes.h"
#include "core/limits.h"
#include "grid/lattice.h"
#include "monitors/spec.h"
#include "geometry/shapes.h"
using fdtd::constants::C0;
using fdtd::constants::INVVACIMP;
using fdtd::constants::I_UNIT;
using fdtd::constants::PI;
using fdtd::constants::PI2;
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;

namespace
{
char _cstring[100]; ///< scratch buffer for the step numbers in file names
const fdtd::monitors::Spec &sp = fdtd::monitors::spec();        ///< what the monitors need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::monitors; // transitional: module state, monitors/state.h

namespace fdtd::monitors
{

// File-local state, formerly globals used only in this file.
namespace
{
FILE *f_qcoeff_sigma = NULL; ///< Output file of the conductivity-based absorption cross section.
FILE *f_qcoeff = NULL;       ///< Output file of the linear Q coefficients.
FILE *f_qcoeff_nl = NULL;    ///< Output file of the harmonic Q coefficients.
/// power through each face [field/harmonic][box: scattering, total][normal axis][near/far][frequency]:
/// this rank's part, and the sum over ranks, root
double q_face_loc[2][2][3][2][fdtd::max_frequencies];
double q_face[2][2][3][2][fdtd::max_frequencies];
double qscattering[fdtd::max_frequencies];    ///< Total scattered power per frequency, outward flux through the scattering box.
double qabsorption[fdtd::max_frequencies];    ///< Total absorbed power per frequency, inward flux through the absorption box.
double qscattering_nl[fdtd::max_frequencies]; ///< Total harmonic scattered power per frequency.
double qabsorption_nl[fdtd::max_frequencies]; ///< Total harmonic absorbed power per frequency.

/**
 * @brief Power through face (n, side) of box b per frequency: the Poynting flux of the DFT fields
 * @verbatim
      P = sum over the face  Re(E_a conj(H_b) - E_b conj(H_a)) d_a d_b,   a = n + 1, b = n + 2, cyclic
   @endverbatim
 * The cells run over the lower tangential axis outside, the higher inside, as before.
 * @return the face area, sum of the cell areas
 */
double face_flux(int set, int b, int n, int side, double Q[])
{
	const auto &F = fdtd::monitors::flux_boxes().box[b].face[n][side];
	const int A = (n + 1) % 3;
	const int B = (n + 2) % 3;
	const int lo = A < B ? A : B;
	const int hi = A < B ? B : A;
	const int range[2][3][2] = {{{subdomain.i_qsca_min, subdomain.i_qsca_max}, {subdomain.j_qsca_min, subdomain.j_qsca_max}, {subdomain.k_qsca_min, subdomain.k_qsca_max}},
								{{subdomain.i_qtot_min, subdomain.i_qtot_max}, {subdomain.j_qtot_min, subdomain.j_qtot_max}, {subdomain.k_qtot_min, subdomain.k_qtot_max}}};
	const double d[3] = {lat.d[0], lat.d[1], lat.d[2]};
	const view4<double> reE_a = F.re[0][set][A];
	const view4<double> imE_a = F.im[0][set][A];
	const view4<double> reE_b = F.re[0][set][B];
	const view4<double> imE_b = F.im[0][set][B];
	const view4<double> reH_a = F.re[1][set][A];
	const view4<double> imH_a = F.im[1][set][A];
	const view4<double> reH_b = F.re[1][set][B];
	const view4<double> imH_b = F.im[1][set][B];
	double q_normalisation_area = 0;
	for (int m = 0; m < sp.nr_frequencies; ++m)
	{
		q_normalisation_area = double(0.);
		double PowerFlux = double(0.);
		for (int p = range[b][lo][0]; p < range[b][lo][1]; ++p)
			for (int q = range[b][hi][0]; q < range[b][hi][1]; ++q)
			{
				q_normalisation_area += d[lo] * d[hi];
				int x[3];
				x[n] = 0, x[lo] = p, x[hi] = q;
				auto at = [&](const view4<double> &v) { return v[m][x[0]][x[1]][x[2]]; };
				const std::complex<double> cEa = at(reE_a) + I_UNIT * at(imE_a);
				const std::complex<double> cEb = at(reE_b) + I_UNIT * at(imE_b);
				const std::complex<double> ccHa = at(reH_a) - I_UNIT * at(imH_a);
				const std::complex<double> ccHb = at(reH_b) - I_UNIT * at(imH_b);
				const std::complex<double> zz = cEa * ccHb - cEb * ccHa;
				PowerFlux += zz.real() * d[lo] * d[hi];
			}
		Q[m] = PowerFlux;
	}
	return q_normalisation_area;
}

bool has_face(int b, int n, int side)
{
	const auto &S = subdomain;
	const bool f[2][3][2] = {{{S.qsca_nea_x, S.qsca_far_x}, {S.qsca_nea_y, S.qsca_far_y}, {S.qsca_nea_z, S.qsca_far_z}},
							 {{S.qtot_nea_x, S.qtot_far_x}, {S.qtot_nea_y, S.qtot_far_y}, {S.qtot_nea_z, S.qtot_far_z}}};
	return f[b][n][side];
}

/// -near_x + far_x - near_y + far_y - near_z + far_z at frequency m, outward flux of the box
double box_flux(const double (&q)[3][2][fdtd::max_frequencies], int m)
{
	return -q[0][0][m] + q[0][1][m] - q[1][0][m] + q[1][1][m] - q[2][0][m] + q[2][1][m];
}
} // namespace

/**
 * @brief Calculates scattering, absorption, and extinction cross-sections according to Mie theory. 
 * This is calculated when TF/SF is enabled. The nanostructure must be enclosed within the TF/SF box. 
 * The method works both for homogeneous and layered background. 
 * For certain geometries, the normalization by the projected area is implemented 
 * e.g., sphere, monopole, dipole, so that the method outputs the scattering,
 * absorption, and extinction efficiencies.
 */
/**
 * @brief Area that normalises the Q coefficients [m^2]
 *
 * The geometric cross section of the structure when its shape declares one, a physics quantity:
 * the projected area from the shape's own dimensions; otherwise the unnormalised flux for the
 * shapes that never had a reference area, rounded triangle, U and H shapes, else the TF/SF box area.
 */
static double compute_q_normalisation_area(double tfsf_area)
{
	if (geometric_cross_section > 0)
		return geometric_cross_section;
	if ((sp.shape == fdtd::geometry::shape::rounded_triangle) || (sp.shape == fdtd::geometry::shape::u_shape) || (sp.shape == fdtd::geometry::shape::h_shape))
		return double(1);
	return tfsf_area;
}

void write_q_coefficients()
{
	int m;
	const int sets = sp.hydrodynamics ? 2 : 1;
	// area_loc: the area of the total box's far y face on this rank, printed as a check
	area_loc = 0;
	for (int set = 0; set < sets; ++set)
		for (int b = 0; b < 2; ++b)
			for (int n = 0; n < 3; ++n)
				for (int side = 0; side < 2; ++side)
					if (has_face(b, n, side))
					{
						const double area = face_flux(set, b, n, side, q_face_loc[set][b][n][side]);
						if (b == 1 && n == 1 && side == 1)
							area_loc = area;
					}

	for (int set = 0; set < sets; ++set)
		for (int b = 0; b < 2; ++b)
			for (int n = 0; n < 3; ++n)
				for (int side = 0; side < 2; ++side)
				{
					MPI_Reduce(q_face_loc[set][b][n][side], q_face[set][b][n][side], sp.nr_frequencies, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
					MPI_Barrier(MPI_COMM_WORLD);
				}

	MPI_Reduce(&area_loc, &q_normalisation_area, 1, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);

	double Area_inc = (lat.inner(2) - 2 * lat.tfsf_inset[2]) * (lat.inner(0) - 2 * lat.tfsf_inset[0]) * lat.d[0] * lat.d[2];

	if (fdtd::grid::world_rank == 0)
		fdtd::io::log("monitors", "Area xz absorption box numeric= %e\t analytic= %e\n", q_normalisation_area, Area_inc);

	if (sp.abs_sigma)
	{
		integrate_absorption(sp.antenna_material);
		MPI_Reduce(abs_sigma_loc, abs_sigma_spectrum, sp.nr_frequencies, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
	}

	if (fdtd::grid::world_rank == ROOTPROCESS)
	{
		std::string strDataPath_Qcoeff = sp.reports_dir + "/cross_sections.txt";
		std::string strDataPath_Qcoeff_NL = sp.reports_dir + "/cross_sections_harmonic.txt";

		f_qcoeff = fdtd::io::open_report(strDataPath_Qcoeff, "");
		if (sp.hydrodynamics)
			f_qcoeff_nl = fdtd::io::open_report(strDataPath_Qcoeff_NL, "");

		if (sp.abs_sigma)
		{
			std::string strDataPath_Qcoeff_sigma = sp.reports_dir + "/cross_sections_sigma.txt";
			f_qcoeff_sigma = fdtd::io::open_report(strDataPath_Qcoeff_sigma, "");
		}
		fprintf(f_qcoeff, "# LAMBDA[nm]\t Qscatt\t Qabs\t Qextin\n");
		if (f_qcoeff_nl)
			fprintf(f_qcoeff_nl, "# LAMBDA[nm]\t Qscatt\t Qabs\t Qextin\n");
		if (f_qcoeff_sigma)
			fprintf(f_qcoeff_sigma, "# LAMBDA[nm]\t Qscatt\t Qabs_sigma\t Qextin\n");

		if (geometric_cross_section > 0)
			fdtd::io::log("monitors", "Q normalised by the geometric cross section %e m^2\n", geometric_cross_section);
		else
			fdtd::io::log("monitors", "Q normalised by %e m^2 (no geometric cross section declared for this shape)\n", compute_q_normalisation_area(Area_inc));
		if (sp.abs_sigma)
			fdtd::io::log("qcoeff", "lambda_nm Qscatt Qabs Qextin Qabs_sigma");
		else
			fdtd::io::log("qcoeff", "lambda_nm Qscatt Qabs Qextin");

		for (m = 0; m < sp.nr_frequencies; ++m)
		{
			qscattering[m] = box_flux(q_face[0][0], m);
			qabsorption[m] = box_flux(q_face[0][1], m);

			if (sp.hydrodynamics)
			{
				qscattering_nl[m] = box_flux(q_face[1][0], m);
				qabsorption_nl[m] = box_flux(q_face[1][1], m);
			}

			if (((sp.thg || sp.shg) && (m < 0.5 * sp.nr_frequencies)) || (!sp.thg && !sp.shg))
			{
				qscattering[m] /= double(fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real() * INVVACIMP); //amplitude_normalization[m]*amplitude_normalization[m]
				if (sp.hydrodynamics)
					qscattering_nl[m] /= double(fdtd::materials::refractive_index(sp.medium_first, 2.0 * fdtd::sources::analysis_omega[m]).real() * INVVACIMP); //amplitude_normalization[m]*amplitude_normalization[m]
				if (sp.draw)
				{
					qscattering[m] /= compute_q_normalisation_area(Area_inc);

					if (sp.hydrodynamics)
					{
						qscattering_nl[m] /= compute_q_normalisation_area(Area_inc);
					}
				}
				qabsorption[m] /= double(fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real() * INVVACIMP); //amplitude_normalization[m]*amplitude_normalization[m]
				if (sp.hydrodynamics)
					qabsorption_nl[m] /= double(fdtd::materials::refractive_index(sp.medium_first, 2.0 * fdtd::sources::analysis_omega[m]).real() * INVVACIMP); //amplitude_normalization[m]*amplitude_normalization[m]
				if (sp.draw)
				{
					qabsorption[m] /= compute_q_normalisation_area(Area_inc);

					if (sp.hydrodynamics)
					{
						qabsorption_nl[m] /= compute_q_normalisation_area(Area_inc);
					}
				}
				if (sp.abs_sigma)
				{
					abs_sigma_spectrum[m] /= double(fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real() * INVVACIMP); //amplitude_normalization[m]*amplitude_normalization[m]
					if (sp.draw)
					{
						abs_sigma_spectrum[m] /= compute_q_normalisation_area(Area_inc);
					}
				}

			}
			if (sp.abs_sigma)
			{
				fprintf(f_qcoeff_sigma, "%e\t%e\t%e\t%e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), qscattering[m], abs_sigma_spectrum[m], qscattering[m] + abs_sigma_spectrum[m]);
				fdtd::io::log("qcoeff", "%e\t%e\t%e\t%e\t%e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), qscattering[m], -qabsorption[m], qscattering[m] - qabsorption[m], abs_sigma_spectrum[m]);
			}
			else
			{
				fdtd::io::log("qcoeff", "%e\t%e\t%e\t%e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), qscattering[m], -qabsorption[m], qscattering[m] - qabsorption[m]);
				if (sp.hydrodynamics)
					fdtd::io::log("qcoeff_nl", "%e\t%e\t%e\t%e\n", (PI2 * C0 / (2.0 * fdtd::sources::analysis_omega[m])) * double(1.e9), qscattering_nl[m], -qabsorption_nl[m], qscattering_nl[m] - qabsorption_nl[m]);
			}
			fprintf(f_qcoeff, "%e\t%e\t%e\t%e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), qscattering[m], -qabsorption[m], qscattering[m] - qabsorption[m]);
			if (sp.hydrodynamics)
				fprintf(f_qcoeff_nl, "%e\t%e\t%e\t%e\n", (PI2 * C0 / (2.0 * fdtd::sources::analysis_omega[m])) * double(1.e9), qscattering_nl[m], -qabsorption_nl[m], qscattering_nl[m] - qabsorption_nl[m]);
		}
		fclose(f_qcoeff);
		if (sp.hydrodynamics)
			fclose(f_qcoeff_nl);

		if (sp.abs_sigma)
			fclose(f_qcoeff_sigma);
	}
	if (fdtd::grid::world_rank == ROOTPROCESS && (sp.thg || sp.shg || sp.hydrodynamics))
	{
		std::string dataPathART = sp.reports_dir + "/harmonic_efficiency.txt";
		f_art = fdtd::io::open_report(dataPathART, "");

		transm_nl_power = double(0.);
		refl_nl_power = double(0.);
		if (sp.nr_frequencies > 2)
		{
			for (m = 0.5 * sp.nr_frequencies; m < sp.nr_frequencies; ++m)
			{
				transm_nl_power += 2 * sp.dt * sp.dt * qscattering[m] * harmonic_step_frequency;//J
				refl_nl_power += 2 * sp.dt * sp.dt * qabsorption[m] * harmonic_step_frequency;//J
			}
			fprintf(f_art, "# Emax[V/m]\t InputFluence[J/m^2]\t InputEnergy[J]\t NLScatEnergy[J]\t NLAbsorpEnergy[J]\t NLScatEff\t NLAbsorpEff\n");
			fprintf(f_art, "%e\t %e\t %e\t %e\t %e\t %e\t %e\n", sp.e_peak, fdtd::sources::input_power, fdtd::sources::input_power * Area_inc, transm_nl_power, refl_nl_power, transm_nl_power / (fdtd::sources::input_power * Area_inc), refl_nl_power / (fdtd::sources::input_power * Area_inc));
		}
		if (sp.nr_frequencies == 2)
		{
			transm_nl_power = qscattering[1];//W
			refl_nl_power = qabsorption[1];//W
			fprintf(f_art, "# Emax[V/m]\t InputIntensity[W/m^2]\t InputPower[W]\t NLScatPower[W]\t NLAbsorpPower[W]\t NLScatEff\t NLAbsorpEff\n");
			fprintf(f_art, "%e\t %e\t %e\t %e\t %e\t %e\t %e\n", sp.e_peak, fdtd::sources::input_power, fdtd::sources::input_power * Area_inc, transm_nl_power, refl_nl_power, transm_nl_power / (fdtd::sources::input_power * Area_inc), refl_nl_power / (fdtd::sources::input_power * Area_inc));
		}
		fclose(f_art);
	}
}
/**
 * @brief THIS FUNCTION MUST BE IN GENERAL AVOIDED. 
 * In the case of TF/SF, reflectance/transmittance cannot be calculated.
 * This function gives an approximation of R and T by calculating the flux of the Poynting vector
 * through only one edge of the six-edges box used for the cross section. R is mimiked by integrating 
 * over the xz-edge of the scattered field box in the backward direction. T is mimiked by integrating
 * over the xz-edge of the total field box in the forward direction.
 * For metasurfaces with periodic boundary conditions enabled, use write_art_spectra(), monitors/art.h.
 * This function may be useful to calculate reflectance and transmittance when simulating a full metasurface.
 */
void write_art_tfsf()
{
	double Area_inc;
	int m;
	Area_inc = (lat.inner(2) - 2 * lat.tfsf_inset[2]) * (lat.inner(0) - 2 * lat.tfsf_inset[0]) * lat.d[0] * lat.d[2];

	// R from the near y face of the scattering box, T from the far y face of the total box
	double area_r_loc = 0;
	double area_r = 0;
	area_loc = 0;
	for (int set = 0; set < (sp.hydrodynamics ? 2 : 1); ++set)
	{
		if (subdomain.qsca_nea_y)
			area_r_loc = face_flux(set, 0, 1, 0, q_face_loc[set][0][1][0]);
		if (subdomain.qtot_far_y)
			area_loc = face_flux(set, 1, 1, 1, q_face_loc[set][1][1][1]);
		if (set == 0)
			// R and T are normalized by the incident power through the face each is integrated over
			MPI_Allreduce(&area_r_loc, &area_r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Reduce(q_face_loc[set][0][1][0], q_face[set][0][1][0], sp.nr_frequencies, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(q_face_loc[set][1][1][1], q_face[set][1][1][1], sp.nr_frequencies, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Allreduce(&area_loc, &q_normalisation_area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if (fdtd::grid::world_rank == 0)
		fdtd::io::log("monitors", "flux faces: reflection %e m^2, transmission %e m^2 (TF box %e m^2)\n", area_r, q_normalisation_area, Area_inc);

	if (sp.abs_sigma)
	{
		integrate_absorption(sp.antenna_material);
		MPI_Reduce(abs_sigma_loc, abs_sigma_spectrum, sp.nr_frequencies, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if (single_photon)
		{
			integrate_absorption(16);
			MPI_Reduce(abs_sigma_loc, abs_sigma_2, sp.nr_frequencies, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	if (fdtd::grid::world_rank == ROOTPROCESS)
	{
		sprintf(_cstring, "%d", int(sp.steps / 10) + count_art * int(sp.steps / 10));
		str_2 = _cstring;
		std::string strDataPath_Qcoeff = sp.reports_dir + "/art_tfsf_" + str_2 + ".txt";
		std::string strDataPath_Qcoeff_NL = sp.reports_dir + "/art_tfsf_harmonic_" + str_2 + ".txt";
		f_qcoeff = fdtd::io::open_report(strDataPath_Qcoeff, "");
		if (sp.hydrodynamics)
			f_qcoeff_nl = fdtd::io::open_report(strDataPath_Qcoeff_NL, "");
		if (sp.abs_sigma)
		{
			std::string strDataPath_Qcoeff_sigma = sp.reports_dir + "/art_tfsf_sigma.txt";
			f_qcoeff_sigma = fdtd::io::open_report(strDataPath_Qcoeff_sigma, "");
		}
		fprintf(f_qcoeff, "# LAMBDA[nm]\t T\t R\t A\n");
		if (sp.hydrodynamics)
			fprintf(f_qcoeff_nl, "# LAMBDA[nm]\t T\t R\t A\n");

		if (sp.abs_sigma)
		{
			if (single_photon)
				fprintf(f_qcoeff_sigma, "# LAMBDA[nm]\t T\t R\t A_met\t A_Si\t Gain\n");
			else
				fprintf(f_qcoeff_sigma, "# LAMBDA[nm]\t T\t R\t A\n");
		}

		for (m = 0; m < sp.nr_frequencies; ++m)
		{
			if (((sp.thg || sp.shg) && (m < 0.5 * sp.nr_frequencies)) || (!sp.thg && !sp.shg))
			{
				qscattering[m] = -q_face[0][0][1][0][m];
				qabsorption[m] = q_face[0][1][1][1][m];

				qscattering[m] /= double(fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real() * INVVACIMP * area_r);
				qabsorption[m] /= double(fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real() * INVVACIMP * q_normalisation_area);

				if (sp.hydrodynamics)
				{
					qscattering_nl[m] = -q_face[1][0][1][0][m];
					qabsorption_nl[m] = q_face[1][1][1][1][m];

					qscattering_nl[m] /= double(fdtd::materials::refractive_index(sp.medium_first, 2.0 * fdtd::sources::analysis_omega[m]).real() * INVVACIMP * area_r);
					qabsorption_nl[m] /= double(fdtd::materials::refractive_index(sp.medium_first, 2.0 * fdtd::sources::analysis_omega[m]).real() * INVVACIMP * q_normalisation_area);
				}

				if (sp.abs_sigma)
				{
					abs_sigma_spectrum[m] /= double(fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real() * INVVACIMP * Area_inc);
					if (single_photon)
						abs_sigma_2[m] /= double(fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real() * INVVACIMP * Area_inc);
				}
			}
			if (sp.abs_sigma)
			{
				if (single_photon)
					fprintf(f_qcoeff_sigma, "%e\t%e\t%e\t%e\t%e\t%e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), qabsorption[m], 1 - qabsorption[m] - abs_sigma_spectrum[m] - abs_sigma_2[m], abs_sigma_spectrum[m], abs_sigma_2[m], abs_sigma_2[m] * Area_inc / (PI * double(0.5E-06) * double(0.5E-06)));
				else
					fprintf(f_qcoeff_sigma, "%e\t%e\t%e\t%e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), qabsorption[m], 1 - qabsorption[m] - abs_sigma_spectrum[m], abs_sigma_spectrum[m]);
			}
			fprintf(f_qcoeff, "%e\t%e\t%e\t%e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), qabsorption[m], qscattering[m], 1 - qscattering[m] - qabsorption[m]);
			if (sp.hydrodynamics)
				fprintf(f_qcoeff_nl, "%e\t%e\t%e\t%e\n", (PI2 * C0 / (2.0 * fdtd::sources::analysis_omega[m])) * double(1.e9), qabsorption_nl[m], qscattering_nl[m], 1 - qscattering_nl[m] - qabsorption_nl[m]);
		}
		fclose(f_qcoeff);
		if (sp.hydrodynamics)
			fclose(f_qcoeff_nl);
		if (sp.abs_sigma)
			fclose(f_qcoeff_sigma);
	}

	if (fdtd::grid::world_rank == ROOTPROCESS && (sp.thg || sp.shg || sp.hydrodynamics))
	{
		std::string dataPathART = sp.reports_dir + "/harmonic_efficiency_art.txt";
		f_art = fdtd::io::open_report(dataPathART, "");

		transm_nl_power = double(0.);
		refl_nl_power = double(0.);
		if (sp.nr_frequencies > 2)
		{
			for (m = 0.5 * sp.nr_frequencies; m < sp.nr_frequencies; ++m)
			{
				transm_nl_power += 2 * sp.dt * sp.dt * (-q_face[0][0][1][0][m]) * harmonic_step_frequency;//J
				refl_nl_power += 2 * sp.dt * sp.dt * q_face[0][1][1][1][m] * harmonic_step_frequency;//J
			}
			fprintf(f_art, "# Emax[V/m]\t InputFluence[J/m^2]\t InputEnergy[J]\t NLScatEnergy[J]\t NLAbsorpEnergy[J]\t NLScatEff\t NLAbsorpEff\n");
			fprintf(f_art, "%e\t %e\t %e\t %e\t %e\t %e\t %e\n", sp.e_peak, fdtd::sources::input_power, fdtd::sources::input_power * Area_inc, transm_nl_power, refl_nl_power, transm_nl_power / (fdtd::sources::input_power * Area_inc), refl_nl_power / (fdtd::sources::input_power * Area_inc));
		}
		if (sp.nr_frequencies == 2)
		{
			transm_nl_power = -q_face[0][0][1][0][1];//W
			refl_nl_power = q_face[0][1][1][1][1];//W
			fprintf(f_art, "# Emax[V/m]\t InputIntensity[W/m^2]\t InputPower[W]\t NLScatPower[W]\t NLAbsorpPower[W]\t NLScatEff\t NLAbsorpEff\n");
			fprintf(f_art, "%e\t %e\t %e\t %e\t %e\t %e\t %e\n", sp.e_peak, fdtd::sources::input_power, fdtd::sources::input_power * Area_inc, transm_nl_power, refl_nl_power, transm_nl_power / (fdtd::sources::input_power * Area_inc), refl_nl_power / (fdtd::sources::input_power * Area_inc));
		}
		fclose(f_art);
	}
}
} // namespace fdtd::monitors
