/**
 * @file monitors/art.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-01-19
 * 
 * @brief File holding the functions for calculating the Transmittance,
 * reflectance and absorbance
 * 
 * @copyright Copyright 2022
 */

#include "monitors/art.h"
#include "io/text.h"
#include "monitors/state.h"
#include "materials/material_map.h"
#include "grid/subdomain.h"
#include "core/views.h"
#include "grid/state.h"
#include "sources/state.h"
#include "materials/materials.h"
#include "monitors/art_planes.h"
#include "monitors/volume_dft.h"
#include "core/limits.h"
#include "grid/lattice.h"
#include "monitors/spec.h"
using fdtd::constants::C0;
using fdtd::constants::EP0;
using fdtd::constants::INVVACIMP;
using fdtd::constants::PI2;
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;
using fdtd::materials::matind_x;
using fdtd::materials::matind_y;
using fdtd::materials::matind_z;

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
FILE *f_art_sigma = NULL;                                   ///< Output file of the conductivity-based ART absorption.
double tran_loc[fdtd::max_frequencies];                          ///< This rank's part of tran.
double refl_loc[fdtd::max_frequencies];                          ///< This rank's part of refl.
double tran[fdtd::max_frequencies];                              ///< Transmitted power per frequency through ART plane 2.
double refl[fdtd::max_frequencies];                              ///< Reflected power per frequency through ART plane 1.
double area_t;                                                ///< Area of the transmission plane in m^2.
double AreaR;                                                 ///< Area of the reflection plane in m^2.
} // namespace

// Absorbance (a), reflectance (r) and trasmittance (t) coefficients for periodic structures, pbc
/**
 * @brief wrapper function for Transmittance and reflectance calculations
 * 
 * @details See Inidividual Transmittance/Reflectance functions for details
 */
void write_art_spectra(int step)
{
	int m;
	double n_factor; //factor that multiplies the transmittance = Re(1/(n+ik)) = n/(n^2+k^2)
	area_loc = 0;
	AreaR = 0;

	// R: the flux back through the near plane; T: the flux through the far plane
	const fdtd::monitors::ArtPlanes &art = fdtd::monitors::art_planes();
	if (subdomain.art_nea_y)
	{
		area_loc = art.flux(0, refl_loc);
		for (m = 0; m < sp.nr_frequencies; ++m)
			refl_loc[m] = -refl_loc[m];
	}

	MPI_Reduce(refl_loc, refl, sp.nr_frequencies, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Reduce(&area_loc, &AreaR, 1, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
	area_loc = 0;
	area_t = 0;

	if (subdomain.art_far_y)
		area_loc = art.flux(1, tran_loc);

	MPI_Reduce(tran_loc, tran, sp.nr_frequencies, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Reduce(&area_loc, &area_t, 1, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);

	if (sp.abs_sigma)
	{
		integrate_absorption(sp.antenna_material);
		MPI_Reduce(abs_sigma_loc, abs_sigma_spectrum, sp.nr_frequencies, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if (single_photon)
		{
			integrate_absorption(sp.medium_second);
			MPI_Reduce(abs_sigma_loc, abs_sigma_2, sp.nr_frequencies, MPI_DOUBLE, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	if (fdtd::grid::world_rank == ROOTPROCESS)
	{

		if (count_art == 9) // the final spectra: step is the last completed step
		{
			sprintf(_cstring, "%d", step);
		}
		else
		{
			sprintf(_cstring, "%d", step);
		}

		str_2 = _cstring;
		std::string dataPathART;
		fdtd::io::log("monitors", "AreaR: %e\t AreaT: %e\n", AreaR, area_t);

		dataPathART = sp.reports_dir + "/ART_" + fdtd::sources::pol + "_" + str_2 + ".txt";

		f_art = fopen(dataPathART.c_str(), "w+");
		fprintf(f_art, "LAMBDA[nm]\t T\t R\t A\n");

		if (sp.abs_sigma)
		{
			std::string strDataPath_ART_sigma = sp.reports_dir + "/ART_" + fdtd::sources::pol + "_sigma.txt";
			f_art_sigma = fopen(strDataPath_ART_sigma.c_str(), "w+");
			if (single_photon)
			{
				fprintf(f_art_sigma, "LAMBDA[nm]\t T\t R\t A_met\t A_Si\t Err\n");
			}
			else
			{
				fprintf(f_art_sigma, "LAMBDA[nm]\t T\t R\t A\t Err\n");
			}
		}

		for (m = 0; m < sp.nr_frequencies; ++m)
		{

			if (((sp.thg || sp.shg) && (m < 0.5 * sp.nr_frequencies)) || (sp.hydrodynamics && (m < (1.0/3.0) * sp.nr_frequencies)) || (!sp.thg && !sp.shg))
			{
				n_factor = fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real();
				//DFT normalization is already done; we now introduce the parameters to finalize the normalization with respect
				//to the energy of the incident pulse
				tran[m] /= (n_factor * INVVACIMP * area_t); //amplitude_normalization[m]*amplitude_normalization[m]);
				refl[m] /= (n_factor * INVVACIMP * AreaR); //amplitude_normalization[m]*amplitude_normalization[m]);

				if (count_art != 9)
				{
					tran[m] /= (fdtd::sources::amplitude_normalization[m] * fdtd::sources::amplitude_normalization[m]);
					refl[m] /= (fdtd::sources::amplitude_normalization[m] * fdtd::sources::amplitude_normalization[m]);
				}
				//media mobile
				//tran[m] = (Tran_prev[m]*count_art+tran[m])/(count_art+1)
				//refl[m] = (Refl_prev[m]*count_art+refl[m])/(count_art+1)
			}
			//		{  //tran and refl are in W. The piece on the right is the power in W of the incident excitation.
			// 			Joule = W*s.
			//			to get energy in the freq. domain, we do dt*dt*tran*d_freq = Joule

			fprintf(f_art, "%e\t %e\t %e\t %e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), tran[m], refl[m], 1 - tran[m] - refl[m]);

			if (sp.abs_sigma)
			{
				if (((sp.thg || sp.shg) && (m < 0.5 * sp.nr_frequencies)) || (!sp.thg && !sp.shg))
				{
					abs_sigma_spectrum[m] /= double(fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real() * INVVACIMP * area_t);
					if (single_photon)
					{
						abs_sigma_2[m] /= double(fdtd::materials::refractive_index(sp.medium_first, fdtd::sources::analysis_omega[m]).real() * INVVACIMP * area_t);
					}
				}
				if (single_photon)
				{
					fprintf(f_art_sigma, "%e\t %e\t %e\t %e\t %e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), tran[m], refl[m], abs_sigma_spectrum[m], abs_sigma_2[m]);
				}
				else
				{
					fprintf(f_art_sigma, "%e\t %e\t %e\t %e\t %e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), tran[m], refl[m], abs_sigma_spectrum[m], 1 - tran[m] - refl[m] - abs_sigma_spectrum[m]);
				}
			}
		}
		fclose(f_art);
		if (sp.abs_sigma)
		{
			fclose(f_art_sigma);
		}
	}

	if (fdtd::grid::world_rank == ROOTPROCESS && (sp.thg || sp.shg || sp.hydrodynamics) && count_art == 9)
	{
		std::string dataPathART = sp.reports_dir + "/NL_Efficiency_PBC.txt";
		f_art = fopen(dataPathART.c_str(), "w+");

		transm_nl_power = double(0.);
		refl_nl_power = double(0.);
		if (sp.nr_frequencies > 2)
		{
			for (m = 0.5 * sp.nr_frequencies; m < sp.nr_frequencies; ++m)
			{
				transm_nl_power += 2 * sp.dt * sp.dt * tran[m] * harmonic_step_frequency;
				refl_nl_power += 2 * sp.dt * sp.dt * refl[m] * harmonic_step_frequency;
			}
			fprintf(f_art, "Emax[V/m]\t InputFluence[J/m^2]\t InputEnergy[J]\t NLScatEnengy[J]\t NLAbsorpEnergy[J]\t NLScatEff\t NLAbsorpEff\n");
			fprintf(f_art, "%e\t %e\t %e\t %e\t %e\t %e\t %e\n", sp.e_peak, fdtd::sources::input_power, fdtd::sources::input_power * AreaR, transm_nl_power, refl_nl_power, transm_nl_power / (fdtd::sources::input_power * AreaR), refl_nl_power / (fdtd::sources::input_power * AreaR));
		}

		if (sp.nr_frequencies == 2)
		{
			transm_nl_power = tran[1];
			refl_nl_power = refl[1];
			fprintf(f_art, "Emax[V/m]\t InputIntensity[W/m^2]\t InputPower[W]\t NLScatPower[W]\t NLAbsorpPower[W]\t NLScatEff\t NLAbsorpEff\n");
			fprintf(f_art, "%e\t %e\t %e\t %e\t %e\t %e\t %e\n", sp.e_peak, fdtd::sources::input_power, fdtd::sources::input_power * AreaR, transm_nl_power, refl_nl_power, transm_nl_power / (fdtd::sources::input_power * AreaR), refl_nl_power / (fdtd::sources::input_power * AreaR));
		}
		fclose(f_art);
	}
}

// T

/**
 * Each cell contributes |E_c|^2 of the components c whose material is `material`, times the cell
 * volume times the number of such components / 3; summed left to right in x, y, z order as before.
 */
void integrate_absorption(int material)
{
	const VolumeDft &V = volume_dft();
	const view4<double> re[3] = {V.re[0][0], V.re[0][1], V.re[0][2]}, im[3] = {V.im[0][0], V.im[0][1], V.im[0][2]};
	const ctensor<int> *map[3] = {&matind_x, &matind_y, &matind_z};
	const double Vol = lat.d[0] * lat.d[1] * lat.d[2];
	for (int m = 0; m < sp.nr_frequencies; ++m)
	{
		double Abs = double(0.);
		for (int i = 0; i < subdomain.inner_dim_x; ++i)
			for (int j = 0; j < subdomain.inner_dim_y; ++j)
				for (int k = 0; k < subdomain.inner_dim_z; ++k)
				{
					double sum = double(0.);
					int count = 0;
					for (int c = 0; c < 3; ++c)
						if (int((*map[c])[i][j][k]) == material)
						{
							sum = count ? sum + re[c][m][i][j][k] * re[c][m][i][j][k] : re[c][m][i][j][k] * re[c][m][i][j][k];
							sum = sum + im[c][m][i][j][k] * im[c][m][i][j][k];
							++count;
						}
					if (count == 3)
						Abs += sum * Vol;
					else if (count == 2)
						Abs += sum * Vol * double(2.) / double(3.);
					else if (count == 1)
						Abs += sum * Vol * double(1.) / double(3.);
				}
		abs_sigma_loc[m] = EP0 * fdtd::sources::analysis_omega[m] * fdtd::materials::complex_permittivity(material, fdtd::sources::analysis_omega[m]).imag() * Abs;
	}
}
} // namespace fdtd::monitors
