/**
 * @file materials/parameters.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-06-15
 * 
 * @brief This file stores the material model parameters.
 * 
 * @details Each material has an index:
 * 0 is vacuum\n
 * 1 - 99 are non-dispersive materials read in from pphmatini.txt\n
 * 
 * Dispersion Offset = 100, 100 is the first dispersive material\n\n
 * 
 * 100 - 149 are Drude + 2CP materials\n
 * 150 - 199 are Drude\n
 * 200 - 249 are Lorentz\n
 * 250 - 299 are Modified Lorentz\n
 * 300 - 349 are Nonlinear Plasmonics, Hydrodynamic PLasma model + 2CP\n
 * 350 - 399 are Drude + Lorentz, with the possibility of nonlinear modelling\n\n
 * 
 * 100 - Gold\n
 * 101 - Silver\n
 * 102 - Graphene\n
 * 103 - Gold-Silver alloy\n
 * 300 - Nonlinear Gold, must be used with Hydrodynamics\n
 * 301 - Nonlinear Silver, must be used with Hydrodynamics\n
 * 302 - Arbitrary Nonlinear Hydroynamics + 2CP materials, must be used with Hydrodynamics\n
 * 304 - Nonlinear Gold, must be used with Hydrodynamics\n
 * 350 - Drude-Lorentz ITO, for multiple thicknesses\n
 * 351 - Arbitrary Drude Lorentz\n
 * 
 * @copyright Copyright 2022
 * 
 */
#include "materials/codes.h"
#include "materials/state.h"
#include "core/constants.h"
#include "materials/nonlinear/state.h"
#include "materials/hydro/state.h"
#include "io/text.h"
#include "materials/thermal/thermal.h"
#include "materials/materials.h"
#include "materials/inputs.h"
#include "grid/lattice.h"
using fdtd::constants::C0;
using fdtd::constants::DEFAULT_TEMPERATURE;
using fdtd::constants::E0;
using fdtd::constants::EP0;
using fdtd::constants::I_UNIT;
using fdtd::constants::ME;
using fdtd::constants::PI;
using fdtd::materials::nr_disp_mat_d;

namespace
{
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace

// File-local state, formerly globals used only in this file.
namespace
{
bool is_chi_2_f_d;                                          ///< chi(2) of the first medium has a nonzero diagonal entry.
bool is_chi_2_f_o;                                          ///< chi(2) of the first medium has a nonzero off-diagonal entry.
bool is_chi_3_f_d;                                          ///< chi(3) of the first medium has a nonzero diagonal entry.
bool is_chi_3_f_o;                                          ///< chi(3) of the first medium has a nonzero off-diagonal entry.
bool is_chi_2_m_d;                                          ///< chi(2) of the middle layer has a nonzero diagonal entry.
bool is_chi_2_m_o;                                          ///< chi(2) of the middle layer has a nonzero off-diagonal entry.
bool is_chi_3_m_d;                                          ///< chi(3) of the middle layer has a nonzero diagonal entry.
bool is_chi_3_m_o;                                          ///< chi(3) of the middle layer has a nonzero off-diagonal entry.
bool is_chi_2_s_d;                                          ///< chi(2) of the second medium has a nonzero diagonal entry.
bool is_chi_2_s_o;                                          ///< chi(2) of the second medium has a nonzero off-diagonal entry.
bool is_chi_3_s_d;                                          ///< chi(3) of the second medium has a nonzero diagonal entry.
bool is_chi_3_s_o;                                          ///< chi(3) of the second medium has a nonzero off-diagonal entry.
double chi_nl[nr_disp_mat_d];                                 ///< Scalar nonlinear coefficient per dispersive material.
double nonlinear_lorentz_beta;                                ///< Nonlinear Lorentz coefficient chi_nl / EP0^2.
} // namespace

/**
 * @brief initializes all relevant global variables for specified materials
 * 
 * @details based off globally available flags, it sets up all required
 * constants from relative permittivities to nonlinear susceptibilities.
 * It does this for the structuring as well as the middle layers. 
 * 
 */
void fdtd::materials::set_parameters(const Inputs &in)
{
	auto &mp = fdtd::materials::model_parameters();
	int m;
	int mm;
	int n;

	//-Izzat-
	//D=Diag, O=Off diagonal
	//A=Antenna, F=First layer
	//M=Middle layer, S=Second layer
	fdtd::materials::nonlinear::is_chi_2_a_d = false;
	fdtd::materials::nonlinear::is_chi_2_a_o = false;
	fdtd::materials::nonlinear::is_chi_3_a_d = false;
	fdtd::materials::nonlinear::is_chi_3_a_o = false;

	is_chi_2_f_d = false;
	is_chi_2_f_o = false;
	is_chi_3_f_d = false;
	is_chi_3_f_o = false;

	is_chi_2_m_d = false;
	is_chi_2_m_o = false;
	is_chi_3_m_d = false;
	is_chi_3_m_o = false;

	is_chi_2_s_d = false;
	is_chi_2_s_o = false;
	is_chi_3_s_d = false;
	is_chi_3_s_o = false;

	for (m = 0; m < NRMAXMATERIALS; ++m)
	{
		mp.epre[m] = double(0.);
		mp.sig_m[m] = double(0.);
		if (m < DispMateOffsetD)
			fdtd::materials::nonlinear::nl_dl[m] = false;
		if (m >= DispMateOffsetD && m < NRMAXMATERIALS)
		{
			mm = m - DispMateOffsetD;
			mp.epre_sta[mm] = double(0.);
			mp.epre_inf[mm] = double(0.);
			mp.ome_plasma[mm] = double(0.);
			mp.gam_plasma[mm] = double(0.);
			chi_nl[mm] = double(0.);
			fdtd::materials::nonlinear::nl[mm] = false;
			for (n = 0; n < NrMaxCP; ++n)
			{
				mp.big_omega_cp[mm][n] = double(0.);
				mp.big_gamma_cp[mm][n] = double(0.);
				mp.a_coeff_cp[mm][n] = double(0.);
				mp.phi_cp[mm][n] = double(0.);
			}
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// dfg_optim: air (0) and GaAs (2) with one index for the optimization in GaAs -> 1
	// graphene_optim: air (0), SiO2(3) and Si (14),
	//   two index for the optimization:  1=SiO2 for the gap region in SiO2 -> first cell in SiO2 at the air/SiO2 interface
	//   								 	2=air for the gap region in air -> last cell in air at the air/SiO2 interface
	//  								 	space not occupied by metal

	//ARIA
	mp.epre[0] = double(1.); //double(1.);
	mp.sig_m[0] = double(0.);

	mp.epre[1] = in.n[1] * in.n[1]; //double(1.7689); //H2O //
	mp.sig_m[1] = double(0.);

	mp.epre[2] = in.n[2] * in.n[2]; //double(12.25); //GaAs
	mp.sig_m[2] = double(0.);

	mp.epre[3] = in.n[3] * in.n[3]; //double(2.3716); //double(2.34); // SiO2;
	mp.sig_m[3] = double(0.);

	for (n = 4; n < 100; n++)
	{
		mp.epre[n] = mp.n_refr[n] * mp.n_refr[n];
		mp.sig_m[n] = double(0.);
	}

	if (in.graphene_optim)
	{
		mp.epre[1] = double(2.3716); //double(2.34); Avg field in SiO2
		mp.sig_m[1] = double(0.);
		mp.epre[2] = double(1); // Avg field in air
		mp.sig_m[2] = double(0.);
		mp.epre[3] = double(2.3716); //SiO2
		mp.sig_m[3] = double(0.);
	}

	//
	// 		NONLINEAR TENSOR IMPLEMENTATION -- IZZAT 2020
	//
	fdtd::materials::nonlinear::is_chi_2_a_d = (in.chi2[0][0] != double(0.) ||
				  in.chi2[0][7] != double(0.) ||
				  in.chi2[0][14] != double(0.)); // at leat one nonzero diagonal element

	fdtd::materials::nonlinear::is_chi_2_a_o = (in.chi2[0][1] == double(0.) &&
				  in.chi2[0][2] == double(0.) &&
				  in.chi2[0][3] == double(0.) &&
				  in.chi2[0][4] == double(0.) &&
				  in.chi2[0][5] == double(0.) &&
				  in.chi2[0][6] == double(0.) &&
				  in.chi2[0][8] == double(0.) &&
				  in.chi2[0][9] == double(0.) &&
				  in.chi2[0][10] == double(0.) &&
				  in.chi2[0][11] == double(0.) &&
				  in.chi2[0][12] == double(0.) &&
				  in.chi2[0][13] == double(0.) &&
				  in.chi2[0][15] == double(0.) &&
				  in.chi2[0][16] == double(0.) &&
				  in.chi2[0][17] == double(0.)); // all offdiagonal elements are zero

	fdtd::materials::nonlinear::is_chi_3_a_d = (in.chi3[0][0] != double(0.) ||
				  in.chi3[0][7] != double(0.) ||
				  in.chi3[0][14] != double(0.)); // at leat one nonzero diagonal element

	fdtd::materials::nonlinear::is_chi_3_a_o = (in.chi3[0][1] == double(0.) &&
				  in.chi3[0][2] == double(0.) &&
				  in.chi3[0][3] == double(0.) &&
				  in.chi3[0][4] == double(0.) &&
				  in.chi3[0][5] == double(0.) &&
				  in.chi3[0][6] == double(0.) &&
				  in.chi3[0][8] == double(0.) &&
				  in.chi3[0][9] == double(0.) &&
				  in.chi3[0][10] == double(0.) &&
				  in.chi3[0][11] == double(0.) &&
				  in.chi3[0][12] == double(0.) &&
				  in.chi3[0][13] == double(0.) &&
				  in.chi3[0][15] == double(0.) &&
				  in.chi3[0][16] == double(0.) &&
				  in.chi3[0][17] == double(0.) &&
				  in.chi3[0][18] == double(0.) &&
				  in.chi3[0][19] == double(0.) &&
				  in.chi3[0][20] == double(0.)); // all offdiagonal elements are zero

	if (in.antenna >= DispMateOffsetD && in.antenna < NRMAXMATERIALS)
	{
		for (m = 0; m < 18; ++m)
		{
			fdtd::materials::nonlinear::chi_2[in.antenna - DispMateOffsetD][m] = in.chi2[0][m]; //chi2 of antenna medium
		}
		for (m = 0; m < 21; ++m)
		{
			fdtd::materials::nonlinear::chi_3[in.antenna - DispMateOffsetD][m] = in.chi3[0][m]; //chi3 of antenna medium
		}
		if (fdtd::materials::nonlinear::is_chi_2_a_d || !fdtd::materials::nonlinear::is_chi_2_a_o || fdtd::materials::nonlinear::is_chi_3_a_d || !fdtd::materials::nonlinear::is_chi_3_a_o)
			fdtd::materials::nonlinear::nl[in.antenna - DispMateOffsetD] = true;
	}
	else if (in.antenna < DispMateOffsetD)
	{
		for (m = 0; m < 18; ++m)
		{
			fdtd::materials::nonlinear::chi2_dl[in.antenna][m] = in.chi2[0][m]; //chi2 of antenna medium
		}
		for (m = 0; m < 21; ++m)
		{
			fdtd::materials::nonlinear::chi3_dl[in.antenna][m] = in.chi3[0][m]; //chi3 of antenna medium
		}
		if (fdtd::materials::nonlinear::is_chi_2_a_d && fdtd::materials::nonlinear::is_chi_2_a_o)
			fdtd::materials::nonlinear::iso_chi_2_dl[in.antenna] = true;
		if (fdtd::materials::nonlinear::is_chi_3_a_d && fdtd::materials::nonlinear::is_chi_3_a_o)
			fdtd::materials::nonlinear::iso_chi_3_dl[in.antenna] = true;
		if (fdtd::materials::nonlinear::is_chi_2_a_d || !fdtd::materials::nonlinear::is_chi_2_a_o || fdtd::materials::nonlinear::is_chi_3_a_d || !fdtd::materials::nonlinear::is_chi_3_a_o)
			fdtd::materials::nonlinear::nl_dl[in.antenna] = true;
	}

	//FIRST_MEDIUM
	is_chi_2_f_d = (in.chi2[1][0] != double(0.) ||
				  in.chi2[1][7] != double(0.) ||
				  in.chi2[1][14] != double(0.)); // at leat one nonzero diagonal element

	is_chi_2_f_o = (in.chi2[1][1] == double(0.) &&
				  in.chi2[1][2] == double(0.) &&
				  in.chi2[1][3] == double(0.) &&
				  in.chi2[1][4] == double(0.) &&
				  in.chi2[1][5] == double(0.) &&
				  in.chi2[1][6] == double(0.) &&
				  in.chi2[1][8] == double(0.) &&
				  in.chi2[1][9] == double(0.) &&
				  in.chi2[1][10] == double(0.) &&
				  in.chi2[1][11] == double(0.) &&
				  in.chi2[1][12] == double(0.) &&
				  in.chi2[1][13] == double(0.) &&
				  in.chi2[1][15] == double(0.) &&
				  in.chi2[1][16] == double(0.) &&
				  in.chi2[1][17] == double(0.)); // all offdiagonal elements are zero

	is_chi_3_f_d = (in.chi3[1][0] != double(0.) ||
				  in.chi3[1][7] != double(0.) ||
				  in.chi3[1][14] != double(0.)); // at leat one nonzero diagonal element

	is_chi_3_f_o = (in.chi3[1][1] == double(0.) &&
				  in.chi3[1][2] == double(0.) &&
				  in.chi3[1][3] == double(0.) &&
				  in.chi3[1][4] == double(0.) &&
				  in.chi3[1][5] == double(0.) &&
				  in.chi3[1][6] == double(0.) &&
				  in.chi3[1][8] == double(0.) &&
				  in.chi3[1][9] == double(0.) &&
				  in.chi3[1][10] == double(0.) &&
				  in.chi3[1][11] == double(0.) &&
				  in.chi3[1][12] == double(0.) &&
				  in.chi3[1][13] == double(0.) &&
				  in.chi3[1][15] == double(0.) &&
				  in.chi3[1][16] == double(0.) &&
				  in.chi3[1][17] == double(0.) &&
				  in.chi3[1][18] == double(0.) &&
				  in.chi3[1][19] == double(0.) &&
				  in.chi3[1][20] == double(0.)); // all offdiagonal elements are zero

	if (in.first >= DispMateOffsetD && in.first < NRMAXMATERIALS)
	{
		for (m = 0; m < 18; ++m)
			fdtd::materials::nonlinear::chi_2[in.first - DispMateOffsetD][m] = in.chi2[1][m]; //chi2 of the first, superstrate, medium

		for (m = 0; m < 21; ++m)
			fdtd::materials::nonlinear::chi_3[in.first - DispMateOffsetD][m] = in.chi3[1][m]; //chi3 of the first, superstrate, medium

		if (is_chi_2_f_d || !is_chi_2_f_o || is_chi_3_f_d || !is_chi_3_f_o)
			fdtd::materials::nonlinear::nl[in.first - DispMateOffsetD] = true;
	}
	else if (in.first < DispMateOffsetD)
	{
		for (m = 0; m < 18; ++m)
			fdtd::materials::nonlinear::chi2_dl[in.first][m] = in.chi2[1][m]; //chi2 of the first, superstrate, medium

		for (m = 0; m < 21; ++m)
			fdtd::materials::nonlinear::chi3_dl[in.first][m] = in.chi3[1][m]; //chi3 of the first, superstrate, medium

		if (is_chi_2_f_d && is_chi_2_f_o)
			fdtd::materials::nonlinear::iso_chi_2_dl[in.first] = true;

		if (is_chi_3_f_d && is_chi_3_f_o)
			fdtd::materials::nonlinear::iso_chi_3_dl[in.first] = true;

		if (is_chi_2_f_d || !is_chi_2_f_o || is_chi_3_f_d || !is_chi_3_f_o)
			fdtd::materials::nonlinear::nl_dl[in.first] = true;
	}

	//MIDDLE_MEDIUM
	is_chi_2_m_d = (in.chi2[2][0] != double(0.) ||
				  in.chi2[2][7] != double(0.) ||
				  in.chi2[2][14] != double(0.)); // at leat one nonzero diagonal element

	is_chi_2_m_o = (in.chi2[2][1] == double(0.) &&
				  in.chi2[2][2] == double(0.) &&
				  in.chi2[2][3] == double(0.) &&
				  in.chi2[2][4] == double(0.) &&
				  in.chi2[2][5] == double(0.) &&
				  in.chi2[2][6] == double(0.) &&
				  in.chi2[2][8] == double(0.) &&
				  in.chi2[2][9] == double(0.) &&
				  in.chi2[2][10] == double(0.) &&
				  in.chi2[2][11] == double(0.) &&
				  in.chi2[2][12] == double(0.) &&
				  in.chi2[2][13] == double(0.) &&
				  in.chi2[2][15] == double(0.) &&
				  in.chi2[2][16] == double(0.) &&
				  in.chi2[2][17] == double(0.)); // all offdiagonal elements are zero

	is_chi_3_m_d = (in.chi3[2][0] != double(0.) ||
				  in.chi3[2][7] != double(0.) ||
				  in.chi3[2][14] != double(0.)); // at leat one nonzero diagonal element

	is_chi_3_m_o = (in.chi3[2][1] == double(0.) &&
				  in.chi3[2][2] == double(0.) &&
				  in.chi3[2][3] == double(0.) &&
				  in.chi3[2][4] == double(0.) &&
				  in.chi3[2][5] == double(0.) &&
				  in.chi3[2][6] == double(0.) &&
				  in.chi3[2][8] == double(0.) &&
				  in.chi3[2][9] == double(0.) &&
				  in.chi3[2][10] == double(0.) &&
				  in.chi3[2][11] == double(0.) &&
				  in.chi3[2][12] == double(0.) &&
				  in.chi3[2][13] == double(0.) &&
				  in.chi3[2][15] == double(0.) &&
				  in.chi3[2][16] == double(0.) &&
				  in.chi3[2][17] == double(0.) &&
				  in.chi3[2][18] == double(0.) &&
				  in.chi3[2][19] == double(0.) &&
				  in.chi3[2][20] == double(0.)); // all offdiagonal elements are zero

	if (in.middle >= DispMateOffsetD && in.middle < NRMAXMATERIALS)
	{
		for (m = 0; m < 18; ++m)
			fdtd::materials::nonlinear::chi_2[in.middle - DispMateOffsetD][m] = in.chi2[2][m]; //chi2 of the middle medium

		for (m = 0; m < 21; ++m)
			fdtd::materials::nonlinear::chi_3[in.middle - DispMateOffsetD][m] = in.chi3[2][m]; //chi3 of the middle medium

		if (is_chi_2_m_d || !is_chi_2_m_o || is_chi_3_m_d || !is_chi_3_m_o)
			fdtd::materials::nonlinear::nl[in.middle - DispMateOffsetD] = true;
	}
	else if (in.middle < DispMateOffsetD)
	{
		for (m = 0; m < 18; ++m)
			fdtd::materials::nonlinear::chi2_dl[in.middle][m] = in.chi2[2][m]; //chi2 of the middle medium

		for (m = 0; m < 21; ++m)
			fdtd::materials::nonlinear::chi3_dl[in.middle][m] = in.chi3[2][m]; //chi3 of the middle medium

		if (is_chi_2_m_d && is_chi_2_m_o)
			fdtd::materials::nonlinear::iso_chi_2_dl[in.middle] = true;

		if (is_chi_3_m_d && is_chi_3_m_o)
			fdtd::materials::nonlinear::iso_chi_3_dl[in.middle] = true;

		if (is_chi_2_m_d || !is_chi_2_m_o || is_chi_3_m_d || !is_chi_3_m_o)
			fdtd::materials::nonlinear::nl_dl[in.middle] = true;
	}

	//SECOND_MEDIUM
	is_chi_2_s_d = (in.chi2[3][0] != double(0.) ||
				  in.chi2[3][7] != double(0.) ||
				  in.chi2[3][14] != double(0.)); // at leat one nonzero diagonal element

	is_chi_2_s_o = (in.chi2[3][1] == double(0.) &&
				  in.chi2[3][2] == double(0.) &&
				  in.chi2[3][3] == double(0.) &&
				  in.chi2[3][4] == double(0.) &&
				  in.chi2[3][5] == double(0.) &&
				  in.chi2[3][6] == double(0.) &&
				  in.chi2[3][8] == double(0.) &&
				  in.chi2[3][9] == double(0.) &&
				  in.chi2[3][10] == double(0.) &&
				  in.chi2[3][11] == double(0.) &&
				  in.chi2[3][12] == double(0.) &&
				  in.chi2[3][13] == double(0.) &&
				  in.chi2[3][15] == double(0.) &&
				  in.chi2[3][16] == double(0.) &&
				  in.chi2[3][17] == double(0.)); // all offdiagonal elements are zero

	is_chi_3_s_d = (in.chi3[3][0] != double(0.) ||
				  in.chi3[3][7] != double(0.) ||
				  in.chi3[3][14] != double(0.)); // at leat one nonzero diagonal element

	is_chi_3_s_o = (in.chi3[3][1] == double(0.) &&
				  in.chi3[3][2] == double(0.) &&
				  in.chi3[3][3] == double(0.) &&
				  in.chi3[3][4] == double(0.) &&
				  in.chi3[3][5] == double(0.) &&
				  in.chi3[3][6] == double(0.) &&
				  in.chi3[3][8] == double(0.) &&
				  in.chi3[3][9] == double(0.) &&
				  in.chi3[3][10] == double(0.) &&
				  in.chi3[3][11] == double(0.) &&
				  in.chi3[3][12] == double(0.) &&
				  in.chi3[3][13] == double(0.) &&
				  in.chi3[3][15] == double(0.) &&
				  in.chi3[3][16] == double(0.) &&
				  in.chi3[3][17] == double(0.) &&
				  in.chi3[3][18] == double(0.) &&
				  in.chi3[3][19] == double(0.) &&
				  in.chi3[3][20] == double(0.)); // all offdiagonal elements are zero

	if (in.second >= DispMateOffsetD && in.second < NRMAXMATERIALS)
	{
		for (m = 0; m < 18; ++m)
			fdtd::materials::nonlinear::chi_2[in.second - DispMateOffsetD][m] = in.chi2[3][m]; //chi2 of the second, substrate, medium

		for (m = 0; m < 21; ++m)
			fdtd::materials::nonlinear::chi_3[in.second - DispMateOffsetD][m] = in.chi3[3][m]; //chi3 of the second, substrate, medium

		if (is_chi_2_s_d || !is_chi_2_s_o || is_chi_3_s_d || !is_chi_3_s_o)
			fdtd::materials::nonlinear::nl[in.second - DispMateOffsetD] = true;
	}
	else if (in.second < DispMateOffsetD)
	{
		for (m = 0; m < 18; ++m)
			fdtd::materials::nonlinear::chi2_dl[in.second][m] = in.chi2[3][m]; //chi2 of the second, substrate, medium

		for (m = 0; m < 21; ++m)
			fdtd::materials::nonlinear::chi3_dl[in.second][m] = in.chi3[3][m]; //chi3 of the second, substrate, medium

		if (is_chi_2_s_d && is_chi_2_s_o)
			fdtd::materials::nonlinear::iso_chi_2_dl[in.second] = true;

		if (is_chi_3_s_d && is_chi_3_s_o)
			fdtd::materials::nonlinear::iso_chi_3_dl[in.second] = true;

		if (is_chi_2_s_d || !is_chi_2_s_o || is_chi_3_s_d || !is_chi_3_s_o)
			fdtd::materials::nonlinear::nl_dl[in.second] = true;
	}

	double FERMI_VELOCITY = double(in.fermi_velocity);
	double D = 3.0;
	double C = 3.0 * D / (D + 2.0);
	double NL_COEFF = sqrt(C / D);

	//
	// 	DRUDE WITH CRITICAL POINTS MATERIALS, Index 0 is material 100
	//

	// Gold Drude + CP (100)
	mp.nr_cp[0] = 2;
	mp.epre_inf[0] = double(1.1431);
	mp.ome_plasma[0] = double(1.3202e16);
	mp.gam_plasma[0] = double(1.0805e14);
	mp.big_omega_cp[0][0] = double(3.8711e15);
	mp.big_omega_cp[0][1] = double(4.1684e15);
	mp.big_gamma_cp[0][0] = double(4.4642e14);
	mp.big_gamma_cp[0][1] = double(2.3555e15);
	mp.a_coeff_cp[0][0] = double(0.26698);
	mp.a_coeff_cp[0][1] = double(3.0834);
	mp.phi_cp[0][0] = double(-1.2371);
	mp.phi_cp[0][1] = double(-1.0968);
	fdtd::materials::hydro::beta_nl[0] = double(1.39E6 * FERMI_VELOCITY * NL_COEFF);
	fdtd::materials::hydro::diff_coeff[0] = double(in.gnor_diffusion);
	fdtd::materials::hydro::equilibrium_electron_density[0] = 5.476415562682574e28;

	// Silver Drude + CP (101)
	mp.nr_cp[1] = 2;
	mp.epre_inf[1] = double(1.4447);
	mp.ome_plasma[1] = double(1.3280e16);
	mp.gam_plasma[1] = double(9.1269e13);

	mp.big_omega_cp[1][0] = double(8.2749e15);
	mp.big_omega_cp[1][1] = double(6.1998e15);
	mp.big_gamma_cp[1][0] = double(5.1770e15);
	mp.big_gamma_cp[1][1] = double(5.4126e14);
	mp.a_coeff_cp[1][0] = double(-1.5951);
	mp.a_coeff_cp[1][1] = double(0.25261);
	mp.phi_cp[1][0] = double(3.1288);
	mp.phi_cp[1][1] = double(-1.5066);
	fdtd::materials::hydro::diff_coeff[1] = double(in.gnor_diffusion);
	fdtd::materials::hydro::beta_nl[1] = double(1.39E6 * FERMI_VELOCITY * NL_COEFF);
	fdtd::materials::hydro::equilibrium_electron_density[1] = 5.541318196761572e28;

	// 300-800 nm

	double graphene_factor = double(0.34E-9) / lat.d[1]; //1; ////mat_z for graphene in flakes
	mp.nr_cp[2] = 4;

	//sostituito in data 31-08-2017 da graphene con eps_inf=1 sempre
	// Drude+4CP
	//
	//
	//
	// 4CP
	//
	//
	//
	mp.epre_inf[2] = double(1); //double(0.147059)*graphene_factor_;
	mp.ome_plasma[2] = double(2137243151133776) * sqrt(in.graphene_stretch * graphene_factor);
	mp.gam_plasma[2] = double(89972692782122.5625);
	mp.big_omega_cp[2][0] = double(8184468641949238);
	mp.big_omega_cp[2][1] = double(6460017572589322);
	mp.big_omega_cp[2][2] = double(5808656148343107);
	mp.big_omega_cp[2][3] = double(776601661856300.75);

	mp.big_gamma_cp[2][0] = double(3362634339149568);
	mp.big_gamma_cp[2][1] = double(769114537286440);
	mp.big_gamma_cp[2][2] = double(2682758040284065);
	mp.big_gamma_cp[2][3] = double(958380663180557.625);

	mp.a_coeff_cp[2][0] = double(5.088825) * in.graphene_stretch * graphene_factor;
	mp.a_coeff_cp[2][1] = double(2.089763) * in.graphene_stretch * graphene_factor;
	mp.a_coeff_cp[2][2] = double(5.765940) * in.graphene_stretch * graphene_factor;
	mp.a_coeff_cp[2][3] = double(19.965192) * in.graphene_stretch * graphene_factor;

	mp.phi_cp[2][0] = double(-2.046461);
	mp.phi_cp[2][1] = double(-0.414203);
	mp.phi_cp[2][2] = double(0.833003);
	mp.phi_cp[2][3] = double(-0.425294);

	//AuAg0.5 added on 04-05-2018, 270-1200nm (103)
	mp.nr_cp[3] = 2;
	mp.epre_inf[3] = double(2.075413);
	mp.ome_plasma[3] = double(13708586572160636);
	mp.gam_plasma[3] = double(254689639317944.22);
	mp.big_omega_cp[3][0] = double(5368271654434062);
	mp.big_omega_cp[3][1] = double(4407683665864418);
	mp.big_gamma_cp[3][0] = double(2611878741320448);
	mp.big_gamma_cp[3][1] = double(827082004377021.75);
	mp.a_coeff_cp[3][0] = double(1.925214);
	mp.a_coeff_cp[3][1] = double(0.243581);
	mp.phi_cp[3][0] = double(-0.584913);
	mp.phi_cp[3][1] = double(-1.971986);

	//sostituito da nuovo graphene con eps_inf=1
	//
	//
	//

	//sostituito da Drude+4CP graphene (300-600nm) on 06-11-2016
	// Alluminio //Alluminio Drude + CP

	//mat = 7; FREE (Drude+CP)
	//AuAg0.5 added on 30-08-2016

	// it became Drude ITO non doped on 30-05-2018, used as undoped ITO in Lidar sims

	//Au, JC [188-1937 nm, 20 steps] (104)
	mp.nr_cp[4] = 2;
	mp.epre_inf[4] = double(1.000001);

	mp.ome_plasma[4] = double(13323154111071286);
	mp.gam_plasma[4] = double(121579950528611.093750);

	mp.big_omega_cp[4][0] = double(4546657998286429);
	mp.big_omega_cp[4][1] = double(4528164722367042);

	mp.big_gamma_cp[4][0] = double(1537670994815246);
	mp.big_gamma_cp[4][1] = double(1225844085119051.25);

	mp.a_coeff_cp[4][0] = double(6.703278);
	mp.a_coeff_cp[4][1] = double(4.243580);

	mp.phi_cp[4][0] = double(-14.162284);
	mp.phi_cp[4][1] = double(-11.329459);

	//Au, JC [188-1937 nm, 20 steps] (105)
	mp.nr_cp[5] = 4;
	mp.epre_inf[5] = double(0.999999);

	mp.ome_plasma[5] = double(14505735125865654.000000);
	mp.gam_plasma[5] = double(313902550046676.562500);

	mp.big_omega_cp[5][0] = double(784444239861788.750000);
	mp.big_omega_cp[5][1] = double(4676557624708066.000000);
	mp.big_omega_cp[5][2] = double(4615297324130980.000000);
	mp.big_omega_cp[5][3] = double(707789567794795.000000);

	mp.big_gamma_cp[5][0] = double(458334898085825.562500);
	mp.big_gamma_cp[5][1] = double(1043515041277192.500000);
	mp.big_gamma_cp[5][2] = double(1380024957985613.000000);
	mp.big_gamma_cp[5][3] = double(8593420695902424.000000);

	mp.a_coeff_cp[5][0] = double(27.427230);
	mp.a_coeff_cp[5][1] = double(3.358332);
	mp.a_coeff_cp[5][2] = double(6.263830);
	mp.a_coeff_cp[5][3] = double(128.189422);

	mp.phi_cp[5][0] = double(-21.997530);
	mp.phi_cp[5][1] = double(-10.687123);
	mp.phi_cp[5][2] = double(-13.701003);
	mp.phi_cp[5][3] = double(-15.647299);

	//Cromium, Sytchkova [262-2500nm, 40 steps] (106)
	mp.nr_cp[6] = 2;
	mp.epre_inf[6] = double(0.999999);

	mp.ome_plasma[6] = double(17803367210498960.000000);
	mp.gam_plasma[6] = double(3772469982659694.000000);

	mp.a_coeff_cp[6][0] = double(4.913459);
	mp.a_coeff_cp[6][1] = double(8.262380);

	mp.big_omega_cp[6][0] = double(3265585771039616.000000);
	mp.big_omega_cp[6][1] = double(1508041794712576.000000);

	mp.big_gamma_cp[6][0] = double(987243800058741.500000);
	mp.big_gamma_cp[6][1] = double(689203619485056.000000);

	mp.phi_cp[6][0] = double(-12.365186);
	mp.phi_cp[6][1] = double(-8.169607);

	//
	// 		DRUDE MATERIALS, START AT 150
	//
	int DrudeShifter = DrudeCP; // This is 50, hence drude materials start 150

	mp.epre_inf[0 + DrudeShifter] = double(4.2345);
	mp.ome_plasma[0 + DrudeShifter] = double(1.6520E15);
	mp.gam_plasma[0 + DrudeShifter] = double(1.7588E14);

	//Drude McPeak gold 1-1.7 um
	mp.epre_inf[0 + DrudeShifter] = double(8.415618);
	mp.ome_plasma[0 + DrudeShifter] = double(14116992901050342);
	mp.gam_plasma[0 + DrudeShifter] = double(48257596409765.203);

	//Drude Yakubovsky gold 0.8-1.8 um

	//mat = 8; Water Drude -> 8-4=4
	//new

	// it become Drude ITO doped on 30-05-2018, until 14-03-2019

	double m_eff = 0.35 * 9.10938356E-31;
	//used until 10-12-2019, simulations for Lidar paper 1
	mp.epre_inf[1 + DrudeShifter] = double(4.2345); //double(3.9);
	mp.ome_plasma[1 + DrudeShifter] = sqrt((in.carrier_density * 1.6021766208E-19 * 1.6021766208E-19) / (EP0 * m_eff));
	mp.gam_plasma[1 + DrudeShifter] = double(1.7588E14); //*1.8E14);

	mp.epre_inf[2 + DrudeShifter] = double(4.2345); // double(3.9);
	mp.ome_plasma[2 + DrudeShifter] = sqrt((in.carrier_density_back * 1.6021766208E-19 * 1.6021766208E-19) / (EP0 * m_eff)); //when Carr_dens_back=3E26 I should get old results
	mp.gam_plasma[2 + DrudeShifter] = double(1.7588E14);												 //*1.8E14);

	//nITO10.5
	// used from 10-12-2019 to simulate experimental ITO by Chaloon for optimization at 900 nm, ITO on Al2O3
	// CARRIER_DENSITY_BACK =  10.543803592417181e+26
	mp.epre_inf[3 + DrudeShifter] = double(3.986568); //double(3.9);
	mp.ome_plasma[3 + DrudeShifter] = sqrt((in.carrier_density * 1.6021766208E-19 * 1.6021766208E-19) / (EP0 * m_eff));
	mp.gam_plasma[3 + DrudeShifter] = double(307237343979061.37); //*1.8E14);

	mp.epre_inf[4 + DrudeShifter] = double(3.986568); // double(3.9);
	mp.ome_plasma[4 + DrudeShifter] = sqrt((in.carrier_density_back * 1.6021766208E-19 * 1.6021766208E-19) / (EP0 * m_eff)); //when Carr_dens_back=3E26 I should get old results
	mp.gam_plasma[4 + DrudeShifter] = double(307237343979061.37);										 //*1.8E14);

	//nITO2.6
	// used from 24-02-2020 to simulate experimental ITO by Chaloon for optimization at 1550 nm, ITO on Al2O3
	// CARRIER_DENSITY_BACK = 2.6216100743790237e+26
	mp.epre_inf[5 + DrudeShifter] = double(3.881872); //double(3.9);
	mp.ome_plasma[5 + DrudeShifter] = sqrt((in.carrier_density * 1.6021766208E-19 * 1.6021766208E-19) / (EP0 * m_eff));
	mp.gam_plasma[5 + DrudeShifter] = double(43609085969994.5); //*1.8E14);

	mp.epre_inf[6 + DrudeShifter] = double(3.881872); // double(3.9);
	mp.ome_plasma[6 + DrudeShifter] = sqrt((in.carrier_density_back * 1.6021766208E-19 * 1.6021766208E-19) / (EP0 * m_eff)); //when Carr_dens_back=3E26 I should get old results
	mp.gam_plasma[6 + DrudeShifter] = double(43609085969994.5);											 //*1.8E14);

	//Drude ITO becomes mat=9 on 06-11-2016
	//sostituito il 03-08-2016 da Drude ITO, usato in Hybrid nanoantenna paper

	//mat=10 becomes Lorentz medium on 20-10-2016, ZnO

	// ZnO Lorentz
	//    EPRE_INF_LOR = 2.797470
	//    EPRE_STAT = 3.686899
	//    OMEGA_RES = 6198529930985674.0
	//    DELTA = 168.0

	//mat =11; FREE, Drude
	// sostituito da Lorentz in 04-08-2016
	//mat = 11; Alumina fitted in 0.2um-1um -> 11-4=7
	//fitted in 350-750 nm, 250 cycle, on 17/06/2017

	int LorentzShifter = Drude;
	mp.epre_sta[0 + LorentzShifter] = double(2.687159);
	mp.epre_inf[0 + LorentzShifter] = double(1.000003);
	mp.ome_plasma[0 + LorentzShifter] = double(19833347701787432);
	mp.gam_plasma[0 + LorentzShifter] = double(1324977.982949);

	// old Alumina
	////	GCOEFF_LORENTZ[8][0] = double(1.);

	//mat = 12; Silicon Lorentz -> 12-4=8
	//sostituito da SiO2 in data 19-11-2015
	////new
	////	GCOEFF_LORENTZ[8][0] = double(1.);

	//SiO2 used for SHG or THG disp
	//mat = 12; SiO2 fitted in 0.2um-1.5um -> 12-4=8
	mp.epre_sta[1 + LorentzShifter] = double(2.093278);
	mp.epre_inf[1 + LorentzShifter] = double(1);
	mp.ome_plasma[1 + LorentzShifter] = double(19958000170839600);
	mp.gam_plasma[1 + LorentzShifter] = double(0.000158);

	// prova chi2 Lorentz

	//sostituito da ZnMgO 0.35 in data 21-04-2017
	//ripristinato in data 02-03-2017
	//sostituito da ZnMgO 0.35 in data 19-10-2016
	//mat = 13; GaAs Lorentz -> 13-4=9

	//ripristinato in data 21-04-2017
	//sostituito da GaAs in data 02-03-2017
	//// mat = 13; ZnMgO Lorentz -> 13-4=9
	mp.epre_sta[2 + LorentzShifter] = double(3.270008);
	mp.epre_inf[2 + LorentzShifter] = double(1.001769);
	//// GCOEFF_LORENTZ[9][0] = double(1.);
	mp.ome_plasma[2 + LorentzShifter] = double(11168632374123130.0);
	mp.gam_plasma[2 + LorentzShifter] = double(1043.075896);

	//// mat = 13; YiG Lorentz -> 203
	mp.epre_sta[3 + LorentzShifter] = double(4.9);
	mp.epre_inf[3 + LorentzShifter] = double(1.0);
	//// GCOEFF_LORENTZ[9][0] = double(1.);
	mp.ome_plasma[3 + LorentzShifter] = double(3.7699112e+15); //600 THz * 2*pi
	mp.gam_plasma[3 + LorentzShifter] = double(0.5*mp.ome_plasma[3 + LorentzShifter]*1e-6);

	int ModifiedLorentzShifter = Lorentz;

	//mat = 14;  Silicon Lorentz modified -> 14-4 = 10
	mp.nr_cp[0 + ModifiedLorentzShifter] = 2;
	mp.epre_inf[0 + ModifiedLorentzShifter] = double(1.);
	mp.big_omega_cp[0 + ModifiedLorentzShifter][0] = double(2 * PI * C0 / (0.292E-6));
	mp.big_omega_cp[0 + ModifiedLorentzShifter][1] = double(2 * PI * C0 / (0.368E-6));
	mp.big_gamma_cp[0 + ModifiedLorentzShifter][0] = double(2 * PI * C0 * 0.425E6);
	mp.big_gamma_cp[0 + ModifiedLorentzShifter][1] = double(2 * PI * C0 * 0.123E6);
	mp.a_coeff_cp[0 + ModifiedLorentzShifter][0] = double(8.93);
	mp.a_coeff_cp[0 + ModifiedLorentzShifter][1] = double(1.855);
	mp.phi_cp[0 + ModifiedLorentzShifter][0] = double(2 * PI * C0 * 0.087E6);
	mp.phi_cp[0 + ModifiedLorentzShifter][1] = double(2 * PI * C0 * 2.678E6);

	//SiO2 by Lorentz mod

	//mat = 15; GaAs Lorentz modified --> 15-4 = 11
	mp.nr_cp[1 + ModifiedLorentzShifter] = 2;
	mp.epre_inf[1 + ModifiedLorentzShifter] = double(1.836773);
	mp.big_omega_cp[1 + ModifiedLorentzShifter][0] = double(7152104570286982.0);
	mp.big_omega_cp[1 + ModifiedLorentzShifter][1] = double(4639193582216234.0);
	mp.big_gamma_cp[1 + ModifiedLorentzShifter][0] = double(773702244543512.0);
	mp.big_gamma_cp[1 + ModifiedLorentzShifter][1] = double(639454757933480.0);
	mp.a_coeff_cp[1 + ModifiedLorentzShifter][0] = double(4.769416);
	mp.a_coeff_cp[1 + ModifiedLorentzShifter][1] = double(4.529038);
	mp.phi_cp[1 + ModifiedLorentzShifter][0] = double(6);
	mp.phi_cp[1 + ModifiedLorentzShifter][1] = double(1962819372613584.0);
	//mat = 16;  Silicon Lorentz modified -> 16-4 = 12
	mp.nr_cp[2 + ModifiedLorentzShifter] = 2;
	mp.epre_inf[2 + ModifiedLorentzShifter] = double(1.);
	mp.big_omega_cp[2 + ModifiedLorentzShifter][0] = double(2 * PI * C0 / (0.292E-6));
	mp.big_omega_cp[2 + ModifiedLorentzShifter][1] = double(2 * PI * C0 / (0.368E-6));
	mp.big_gamma_cp[2 + ModifiedLorentzShifter][0] = double(2 * PI * C0 * 0.425E6);
	mp.big_gamma_cp[2 + ModifiedLorentzShifter][1] = double(2 * PI * C0 * 0.123E6);
	mp.a_coeff_cp[2 + ModifiedLorentzShifter][0] = double(8.93);
	mp.a_coeff_cp[2 + ModifiedLorentzShifter][1] = double(1.855);
	mp.phi_cp[2 + ModifiedLorentzShifter][0] = double(2 * PI * C0 * 0.087E6);
	mp.phi_cp[2 + ModifiedLorentzShifter][1] = double(2 * PI * C0 * 2.678E6);

	int NonLinFreeElectronShifter = ModLorentz;
	// mat = 17 --> 17-4=13; Gold Drude + CP
	// Oro/Gold

	mp.nr_cp[0 + NonLinFreeElectronShifter] = 2;
	mp.epre_inf[0 + NonLinFreeElectronShifter] = double(1.1431);
	mp.ome_plasma[0 + NonLinFreeElectronShifter] = double(1.3202e16);
	mp.gam_plasma[0 + NonLinFreeElectronShifter] = double(1.0805e14);
	mp.big_omega_cp[0 + NonLinFreeElectronShifter][0] = double(3.8711e15);
	mp.big_omega_cp[0 + NonLinFreeElectronShifter][1] = double(4.1684e15);
	mp.big_gamma_cp[0 + NonLinFreeElectronShifter][0] = double(4.4642e14);
	mp.big_gamma_cp[0 + NonLinFreeElectronShifter][1] = double(2.3555e15);
	mp.a_coeff_cp[0 + NonLinFreeElectronShifter][0] = double(0.26698);
	mp.a_coeff_cp[0 + NonLinFreeElectronShifter][1] = double(3.0834);
	mp.phi_cp[0 + NonLinFreeElectronShifter][0] = double(-1.2371);
	mp.phi_cp[0 + NonLinFreeElectronShifter][1] = double(-1.0968);
	fdtd::materials::hydro::beta_nl[0 + NonLinFreeElectronShifter] = double(FERMI_VELOCITY * NL_COEFF);
	fdtd::materials::hydro::equilibrium_electron_density[0 + NonLinFreeElectronShifter] = 5.476415562682574e28;
	mp.epre[0 + NonLinFreeElectronShifter] = 1.0;
	fdtd::materials::hydro::diff_coeff[0 + NonLinFreeElectronShifter] = double(in.gnor_diffusion);
	fdtd::materials::hydro::effective_electron_mass[0 + NonLinFreeElectronShifter] = E0 * E0 * fdtd::materials::hydro::equilibrium_electron_density[0 + NonLinFreeElectronShifter] / (mp.ome_plasma[0 + NonLinFreeElectronShifter] * mp.ome_plasma[0 + NonLinFreeElectronShifter] * EP0);
	fdtd::materials::hydro::effective_electron_mass[0 + NonLinFreeElectronShifter] = ME;
	//mat = 18 --> 18-4=14; Silver Drude + CP

	// Argento/Silver

	mp.nr_cp[1 + NonLinFreeElectronShifter] = 2;
	mp.epre_inf[1 + NonLinFreeElectronShifter] = double(1.4447);
	mp.ome_plasma[1 + NonLinFreeElectronShifter] = double(1.3280e16);
	mp.gam_plasma[1 + NonLinFreeElectronShifter] = double(9.1269e13);
	mp.big_omega_cp[1 + NonLinFreeElectronShifter][0] = double(8.2749e15);
	mp.big_omega_cp[1 + NonLinFreeElectronShifter][1] = double(6.1998e15);
	mp.big_gamma_cp[1 + NonLinFreeElectronShifter][0] = double(5.1770e15);
	mp.big_gamma_cp[1 + NonLinFreeElectronShifter][1] = double(5.4126e14);
	mp.a_coeff_cp[1 + NonLinFreeElectronShifter][0] = double(-1.5951);
	mp.a_coeff_cp[1 + NonLinFreeElectronShifter][1] = double(0.25261);
	mp.phi_cp[1 + NonLinFreeElectronShifter][0] = double(3.1288);
	mp.phi_cp[1 + NonLinFreeElectronShifter][1] = double(-1.5066);
	fdtd::materials::hydro::beta_nl[1 + NonLinFreeElectronShifter] = double(FERMI_VELOCITY * NL_COEFF);
	fdtd::materials::hydro::equilibrium_electron_density[1 + NonLinFreeElectronShifter] = 5.541318196761572e28;
	mp.epre[1 + NonLinFreeElectronShifter] = 1.0;
	fdtd::materials::hydro::effective_electron_mass[1 + NonLinFreeElectronShifter] = E0 * E0 * fdtd::materials::hydro::equilibrium_electron_density[1 + NonLinFreeElectronShifter] / (mp.ome_plasma[1 + NonLinFreeElectronShifter] * mp.ome_plasma[1 + NonLinFreeElectronShifter] * EP0);

	// mat = 19 --> 19-4=15; Arbitrary Nonlinear Metal Electron

	mp.nr_cp[2 + NonLinFreeElectronShifter] = 2;
	mp.epre_inf[2 + NonLinFreeElectronShifter] = in.custom.eps_inf;
	mp.ome_plasma[2 + NonLinFreeElectronShifter] = in.custom.plasma_frequency;
	mp.gam_plasma[2 + NonLinFreeElectronShifter] = in.custom.plasma_damping;
	mp.big_omega_cp[2 + NonLinFreeElectronShifter][0] = in.custom.cp_frequency;
	mp.big_omega_cp[2 + NonLinFreeElectronShifter][1] = in.custom.cp_frequency_2;
	mp.big_gamma_cp[2 + NonLinFreeElectronShifter][0] = in.custom.cp_damping;
	mp.big_gamma_cp[2 + NonLinFreeElectronShifter][1] = in.custom.cp_damping_2;
	mp.a_coeff_cp[2 + NonLinFreeElectronShifter][0] = in.custom.cp_amplitude;
	mp.a_coeff_cp[2 + NonLinFreeElectronShifter][1] = in.custom.cp_amplitude_2;
	mp.phi_cp[2 + NonLinFreeElectronShifter][0] = in.custom.cp_phase;
	mp.phi_cp[2 + NonLinFreeElectronShifter][1] = in.custom.cp_phase_2;
	fdtd::materials::hydro::beta_nl[2 + NonLinFreeElectronShifter] = double(in.custom.fermi_velocity * FERMI_VELOCITY * NL_COEFF);
	fdtd::materials::hydro::equilibrium_electron_density[2 + NonLinFreeElectronShifter] = in.custom.electron_density;
	fdtd::materials::hydro::effective_electron_mass[2 + NonLinFreeElectronShifter] = E0 * E0 * fdtd::materials::hydro::equilibrium_electron_density[2 + NonLinFreeElectronShifter] / (mp.ome_plasma[2 + NonLinFreeElectronShifter] * mp.ome_plasma[2 + NonLinFreeElectronShifter] * EP0);

	//Au, JC, 20 steps//Izzat//mat index 304
	//188-1937 nm

	mp.nr_cp[4 + NonLinFreeElectronShifter] = 2;
	mp.epre_inf[4 + NonLinFreeElectronShifter] = double(1.000001);
	mp.ome_plasma[4 + NonLinFreeElectronShifter] = double(13323154111071286);
	mp.gam_plasma[4 + NonLinFreeElectronShifter] = double(121579950528611.093750);
	mp.big_omega_cp[4 + NonLinFreeElectronShifter][0] = double(4546657998286429);
	mp.big_omega_cp[4 + NonLinFreeElectronShifter][1] = double(4528164722367042);
	mp.big_gamma_cp[4 + NonLinFreeElectronShifter][0] = double(1537670994815246);
	mp.big_gamma_cp[4 + NonLinFreeElectronShifter][1] = double(1225844085119051.25);
	mp.a_coeff_cp[4 + NonLinFreeElectronShifter][0] = double(6.703278);
	mp.a_coeff_cp[4 + NonLinFreeElectronShifter][1] = double(4.243580);
	mp.phi_cp[4 + NonLinFreeElectronShifter][0] = double(-14.162284);
	mp.phi_cp[4 + NonLinFreeElectronShifter][1] = double(-11.329459);
	fdtd::materials::hydro::beta_nl[4 + NonLinFreeElectronShifter] = double(FERMI_VELOCITY * NL_COEFF);
	fdtd::materials::hydro::equilibrium_electron_density[4 + NonLinFreeElectronShifter] = 5.577390364720603e+28;
	mp.epre[4 + NonLinFreeElectronShifter] = 1.0;
	fdtd::materials::hydro::diff_coeff[4 + NonLinFreeElectronShifter] = double(in.gnor_diffusion);
	fdtd::materials::hydro::effective_electron_mass[4 + NonLinFreeElectronShifter] = ME;

	// mat = 19 --> 19-4=15; Gold Free Electron
	// Oro/Gold

	int DrudeLorentzShifter = NonLinFreeElectron;
	//mat = 21 --> 21-4=17 Nonlinear ITO

	if (in.middle_layer < 50e-9)
	{
		mp.nr_cp[DrudeLorentzShifter + 0] = 1;
		mp.ome_plasma[DrudeLorentzShifter + 0] = 2.6594e15;
		mp.ome_plasma[DrudeLorentzShifter + 0] = 2.63e15;
		mp.ome_plasma[DrudeLorentzShifter + 0] = 2.61e15;

		mp.epre_inf[DrudeLorentzShifter + 0] = double(3.1178);

		mp.plasma_freq_add_term = mp.ome_plasma[17] - fdtd::materials::thermal::plasma_frequency(DEFAULT_TEMPERATURE);
		mp.gam_plasma[DrudeLorentzShifter + 0] = double(2.31806e14);
		mp.big_omega_cp[DrudeLorentzShifter + 0][0] = double(7.42969e15);
		mp.big_gamma_cp[DrudeLorentzShifter + 0][0] = double(5.34217e14);
		mp.a_coeff_cp[DrudeLorentzShifter + 0][0] = double(mp.big_omega_cp[DrudeLorentzShifter + 0][0] * sqrt(3.8466980 - mp.epre_inf[DrudeLorentzShifter + 0]));
		mp.epre[DrudeLorentzShifter + 0] = 1.0;
		chi_nl[DrudeLorentzShifter + 0] = 2.16e-18;
		fdtd::materials::nonlinear::nl[DrudeLorentzShifter + 0] = false;
		nonlinear_lorentz_beta = chi_nl[DrudeLorentzShifter + 0] / pow(EP0, 2);
	}

	else if (in.middle_layer < 300e-9)
	{
		mp.nr_cp[DrudeLorentzShifter + 0] = 1;
		mp.epre_inf[DrudeLorentzShifter + 0] = double(3.404);
		mp.ome_plasma[DrudeLorentzShifter + 0] = 2.8790e+15;
		mp.plasma_freq_add_term = mp.ome_plasma[DrudeLorentzShifter + 0] - fdtd::materials::thermal::plasma_frequency(DEFAULT_TEMPERATURE);
		mp.gam_plasma[DrudeLorentzShifter + 0] = double(1.6712e+14);
		mp.big_omega_cp[DrudeLorentzShifter + 0][0] = double(0.0);
		mp.big_gamma_cp[DrudeLorentzShifter + 0][0] = double(0.0);
		mp.a_coeff_cp[DrudeLorentzShifter + 0][0] = double(0.0);
		mp.epre[DrudeLorentzShifter + 0] = 1.0;
		chi_nl[DrudeLorentzShifter + 0] = 2.16e-18;
		fdtd::materials::nonlinear::nl[DrudeLorentzShifter + 0] = false;
		nonlinear_lorentz_beta = chi_nl[DrudeLorentzShifter + 0] / pow(EP0, 2);
	}
	else
	{
		mp.nr_cp[DrudeLorentzShifter + 0] = 1;
		mp.epre_inf[DrudeLorentzShifter + 0] = double(3.8055);
		mp.ome_plasma[DrudeLorentzShifter + 0] = 2.97e+15;

		mp.plasma_freq_add_term = mp.ome_plasma[DrudeLorentzShifter + 0] - fdtd::materials::thermal::plasma_frequency(DEFAULT_TEMPERATURE);

		mp.gam_plasma[DrudeLorentzShifter + 0] = double(1.3909e+14);

		mp.big_omega_cp[DrudeLorentzShifter + 0][0] = double(0.0);
		mp.big_gamma_cp[DrudeLorentzShifter + 0][0] = double(0.0);
		mp.a_coeff_cp[DrudeLorentzShifter + 0][0] = double(0.0);
		mp.epre[DrudeLorentzShifter + 0] = 1.0;
		chi_nl[DrudeLorentzShifter + 0] = 2.16e-18;
		fdtd::materials::nonlinear::nl[DrudeLorentzShifter + 0] = false;
		nonlinear_lorentz_beta = chi_nl[DrudeLorentzShifter + 0] / pow(EP0, 2);

	}
	fdtd::materials::hydro::equilibrium_electron_density[DrudeLorentzShifter + 0] = 1.5477984226363534e+27;
	//mat = 22 --> 22-4=18 Custum Drude Lorentz Material
	mp.nr_cp[DrudeLorentzShifter + 1] = 1;
	mp.epre_inf[DrudeLorentzShifter + 1] = in.custom.eps_inf;
	mp.ome_plasma[DrudeLorentzShifter + 1] = in.custom.plasma_frequency;
	mp.gam_plasma[DrudeLorentzShifter + 1] = in.custom.plasma_damping;
	mp.big_omega_cp[DrudeLorentzShifter + 1][0] = in.custom.cp_frequency;
	mp.big_gamma_cp[DrudeLorentzShifter + 1][0] = in.custom.cp_damping;
	mp.a_coeff_cp[DrudeLorentzShifter + 1][0] = in.custom.cp_amplitude;
	mp.epre[DrudeLorentzShifter + 1] = 1.0;
	fdtd::materials::nonlinear::nl[DrudeLorentzShifter + 1] = false;
	chi_nl[DrudeLorentzShifter + 1] = 2.16e-18;

	int MagneticLorentzShifter = NonlinearDrudeLorentz;

	//// YiG Lorentz -> 400
	mp.epre_sta[0 + MagneticLorentzShifter] = double(4.9);
	mp.epre_inf[0 + MagneticLorentzShifter] = double(1.0);
	mp.ome_plasma[0 + MagneticLorentzShifter] = double(3.7699112e+15); //600 THz * 2*pi
	mp.gam_plasma[0 + MagneticLorentzShifter] = double(mp.ome_plasma[0 + MagneticLorentzShifter]*1e-6);
	mp.a_3[0 + MagneticLorentzShifter] = double(-2.25e22);

}

/**
 * @brief calculates the fractive index of a material at
 * a given value of omega.
 *
 * @param mat index of the material being evaluated, gives model info
 * @param omega angular frequency to evaluate at
 *
 * @return std::complex<double> refractive index at the given value of omega.
 *
 */
std::complex<double> fdtd::materials::refractive_index(int mat, double omega)
{
	const auto &mp = fdtd::materials::model_parameters();
	int mmm;
	std::complex<double> eps = std::complex<double>(0.);
	std::complex<double> refr = std::complex<double>(0.);
	if (mat >= DispMateOffsetD)
	{
		mmm = mat - DispMateOffsetD;

		if (mmm < DrudeCP)
		{
			//DRUDE + CP
			eps = mp.epre_inf[mmm] - mp.ome_plasma[mmm] * mp.ome_plasma[mmm] / omega / (omega + I_UNIT * mp.gam_plasma[mmm]);
			for (int s = 0; s < mp.nr_cp[mmm]; ++s)
			{

				eps += mp.a_coeff_cp[mmm][s] * mp.big_omega_cp[mmm][s] * (std::complex<double>(exp(std::complex<double>(+I_UNIT * mp.phi_cp[mmm][s]))) / (mp.big_omega_cp[mmm][s] - omega - I_UNIT * mp.big_gamma_cp[mmm][s]) + std::complex<double>(exp(std::complex<double>(-I_UNIT * mp.phi_cp[mmm][s]))) / (mp.big_omega_cp[mmm][s] + omega + I_UNIT * mp.big_gamma_cp[mmm][s]));
			}
		}
		else if (mmm >= DrudeCP && mmm < Drude)
		{ //DRUDE

			eps = mp.epre_inf[mmm] - mp.ome_plasma[mmm] * mp.ome_plasma[mmm] / omega / (omega + I_UNIT * mp.gam_plasma[mmm]);
		}
		else if (mmm >= Drude && mmm < Lorentz)
		{
			//LORENTZ

			eps = mp.ome_plasma[mmm] * mp.ome_plasma[mmm] / (mp.ome_plasma[mmm] * mp.ome_plasma[mmm] - double(2.) * I_UNIT * omega * mp.gam_plasma[mmm] - omega * omega);
			eps *= (mp.epre_sta[mmm] - mp.epre_inf[mmm]);
			eps += mp.epre_inf[mmm];
		}
		else if (mmm >= Lorentz && mmm < ModLorentz)
		{
			//LORENTZ MODIFIED
			eps = mp.epre_inf[mmm];
			for (int s = 0; s < mp.nr_cp[mmm]; ++s)
			{
				eps += mp.a_coeff_cp[mmm][s] * (mp.big_omega_cp[mmm][s] * mp.big_omega_cp[mmm][s] - I_UNIT * omega * mp.phi_cp[mmm][s]) / (mp.big_omega_cp[mmm][s] * mp.big_omega_cp[mmm][s] - double(2.) * I_UNIT * omega * mp.big_gamma_cp[mmm][s] - omega * omega);
			}
		}
	}
	else
	{
		eps = std::complex<double>(mp.epre[mat], 0.);
	}
	refr = std::complex<double>(sqrt((sqrt(eps.real() * eps.real() + eps.imag() * eps.imag()) + eps.real()) / 2), sqrt((sqrt(eps.real() * eps.real() + eps.imag() * eps.imag()) - eps.real()) / 2));

	return refr;
}

/**
 * @brief calculates the permattivity of a material at
 * a given value of omega.
 *
 * @param mat index of the material being evaluated, gives model info
 * @param omega angular frequency to evaluate at
 *
 * @return std::complex<double> permittivity at the given value of omega.
 *
 */
std::complex<double> fdtd::materials::complex_permittivity(int mat, double omega)
{
	const auto &mp = fdtd::materials::model_parameters();
	int mmm;
	std::complex<double> eps = std::complex<double>(0.);
	if (mat >= DispMateOffsetD)
	{
		mmm = mat - DispMateOffsetD;

		if (mmm < DrudeCP)
		{
			//DRUDE + CP
			eps = mp.epre_inf[mmm] - mp.ome_plasma[mmm] * mp.ome_plasma[mmm] / omega / (omega + I_UNIT * mp.gam_plasma[mmm]);
			for (int s = 0; s < mp.nr_cp[mmm]; ++s)
			{

				eps += mp.a_coeff_cp[mmm][s] * mp.big_omega_cp[mmm][s] * (std::complex<double>(exp(std::complex<double>(+I_UNIT * mp.phi_cp[mmm][s]))) / (mp.big_omega_cp[mmm][s] - omega - I_UNIT * mp.big_gamma_cp[mmm][s]) + std::complex<double>(exp(std::complex<double>(-I_UNIT * mp.phi_cp[mmm][s]))) / (mp.big_omega_cp[mmm][s] + omega + I_UNIT * mp.big_gamma_cp[mmm][s]));
			}
		}
		else if (mmm >= DrudeCP && mmm < Drude)
		{ //DRUDE

			eps = mp.epre_inf[mmm] - mp.ome_plasma[mmm] * mp.ome_plasma[mmm] / omega / (omega + I_UNIT * mp.gam_plasma[mmm]);
		}
		else if (mmm >= Drude && mmm < Lorentz)
		{
			//LORENTZ

			eps = mp.ome_plasma[mmm] * mp.ome_plasma[mmm] / (mp.ome_plasma[mmm] * mp.ome_plasma[mmm] - I_UNIT * omega * mp.gam_plasma[mmm] - omega * omega);
			eps *= (mp.epre_sta[mmm] - mp.epre_inf[mmm]);
			eps += mp.epre_inf[mmm];
		}
		else if (mmm >= Lorentz && mmm < ModLorentz)
		{
			//LORENTZ MODIFIED
			eps = mp.epre_inf[mmm];
			for (int s = 0; s < mp.nr_cp[mmm]; ++s)
			{
				eps += mp.a_coeff_cp[mmm][s] * (mp.big_omega_cp[mmm][s] * mp.big_omega_cp[mmm][s] - I_UNIT * omega * mp.phi_cp[mmm][s]) / (mp.big_omega_cp[mmm][s] * mp.big_omega_cp[mmm][s] - double(2.) * I_UNIT * omega * mp.big_gamma_cp[mmm][s] - omega * omega);
			}
		}
	}
	else
	{
		eps = std::complex<double>(mp.epre[mat], 0.);
	}

	return eps;
}

