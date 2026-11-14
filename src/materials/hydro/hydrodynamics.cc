/**
 * @file materials/hydro/hydrodynamics.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-03-01
 *
 * @brief Contains functions required for hydrodynamic plasma model
 *
 * @details This file is the contrains functions for updating the current
 * densities and electron densities. It also contains the functions responsible
 * for populating the EPSINF arrays, required in fdtd::materials::hydro::update_electron_density()
 *
 * @copyright Copyright 2022
 */

#include "io/text.h"
#include "materials/models.h"
#include "materials/hydro/state.h"
#include "materials/thermal/thermal.h"
#include "materials/hydro/hydro.h"
#include "grid/halo_exchange.h"
#include "grid/subdomain.h"
#include "grid/lattice.h"
#include "materials/spec.h"
#include "grid/stretching.h"
using fdtd::constants::E0;
using fdtd::constants::EP0;
using fdtd::constants::MU0;
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::fields::Hx;
using fdtd::fields::Hy;
using fdtd::fields::Hz;
using fdtd::fields::ex_nl;
using fdtd::fields::ey_nl;
using fdtd::fields::ez_nl;
using fdtd::fields::hx_prev;
using fdtd::fields::hy_prev;
using fdtd::fields::hz_prev;
using fdtd::grid::subdomain;
using fdtd::materials::DispMateOffsetD;
using fdtd::materials::ModLorentz;
using fdtd::materials::NRMAXMATERIALS;
using fdtd::materials::NonLinFreeElectron;
using fdtd::materials::NonlinearDrudeLorentz;
using fdtd::materials::matind_x;
using fdtd::materials::matind_y;
using fdtd::materials::matind_z;
using fdtd::materials::background_y;

namespace
{
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::materials::hydro; // transitional: module state, hydro/state.h

namespace
{
} // namespace

/**
 * @brief updates the electron density fields nd_x for hydrodynamic modelling
 *
 * @details This function starts by ensuring that all processes have the
 * nessecary field data, via the fdtd::grid::exchange_margins functions. It does this by using
 * additional field arrays like ex_2 and px_cp_2. The electron density is calculated
 * using Gauss's law. It then finally uses Multiswitch to ensure all processes
 * have the correct values of the electron density at the process boundaries
 * using fdtd::grid::exchange_margins
 *
 */

void fdtd::materials::hydro::update_electron_density() {
  int i, j, k, ii, jj, kk;
  for (i = 0; i <= subdomain.cells_x; i++) {
    ii = i + 1;
    for (j = 0; j <= subdomain.cells_y; j++) {
      jj = j + 1;
      for (k = 0; k <= subdomain.cells_z; k++) {
        kk = k + 1;
        ex_2[ii][jj][kk] = Ex[i][j][k];
        ey_2[ii][jj][kk] = Ey[i][j][k];
        ez_2[ii][jj][kk] = Ez[i][j][k];
        if (fdtd::materials::model_spec().harmonic) {
          ex_2_nl[ii][jj][kk] = ex_nl[i][j][k];
          ey_2_nl[ii][jj][kk] = ey_nl[i][j][k];
          ez_2_nl[ii][jj][kk] = ez_nl[i][j][k];
        }
      }
    }
  }

	fdtd::grid::exchange_margins(ex_2, ey_2, ez_2);
	fdtd::grid::exchange_margins(px_cp_2, py_cp_2, pz_cp_2);
	if (fdtd::materials::model_spec().harmonic)
	{
		fdtd::grid::exchange_margins(ex_2_nl, ey_2_nl, ez_2_nl);
		fdtd::grid::exchange_margins(px_cp_2_nl, py_cp_2_nl, pz_cp_2_nl);
	}

  if (!subdomain.cpml_only) 
  {
	for (i = 0; i <= subdomain.cells_x; i++) 
	{
      ii = i + 1;
      for (j = 0; j <= subdomain.cells_y; j++) 
	  {
        jj = j + 1;
        for (k = 0; k <= subdomain.cells_z; k++) 
		{
          kk = k + 1;
          nd_x[ii][jj][kk] = (-1.0 / E0) * (EP0 * (lat.inv_d[0] * (ex_2[ii][jj][kk] * epsinf_x[ii][jj][kk] - ex_2[ii - 1][jj][kk] * epsinf_x[ii - 1][jj][kk]) +
              lat.inv_d[1] * (ey_2[ii][jj][kk] * epsinf_y[ii][jj][kk] - ey_2[ii][jj - 1][kk] * epsinf_y[ii][jj - 1][kk]) +
              lat.inv_d[2] * (ez_2[ii][jj][kk] * epsinf_z[ii][jj][kk] - ez_2[ii][jj][kk - 1] * epsinf_z[ii][jj][kk - 1])) +
            lat.inv_d[0] * (px_cp_2[ii][jj][kk] - px_cp_2[ii - 1][jj][kk]) +
            lat.inv_d[1] * (py_cp_2[ii][jj][kk] - py_cp_2[ii][jj - 1][kk]) +
            lat.inv_d[2] * (pz_cp_2[ii][jj][kk] - pz_cp_2[ii][jj][kk - 1]));

					nd_z[ii][jj][kk] = (-1.0 / E0) * (EP0 * (lat.inv_d[0] * (ex_2[ii][jj][kk] - ex_2[ii - 1][jj][kk]) + lat.inv_d[1] * (ey_2[ii][jj][kk] - ey_2[ii][jj - 1][kk]) + lat.inv_d[2] * (ez_2[ii][jj][kk] - ez_2[ii][jj][kk - 1])));

					if (fdtd::materials::model_spec().harmonic)
					{
						nd_y[ii][jj][kk] = (-1.0 / E0) * (EP0 * (lat.inv_d[0] * (ex_2_nl[ii][jj][kk] * epsinf_x[ii][jj][kk] - ex_2_nl[ii - 1][jj][kk] * epsinf_x[ii - 1][jj][kk]) +
																lat.inv_d[1] * (ey_2_nl[ii][jj][kk] * epsinf_y[ii][jj][kk] - ey_2_nl[ii][jj - 1][kk] * epsinf_y[ii][jj - 1][kk]) +
																lat.inv_d[2] * (ez_2_nl[ii][jj][kk] * epsinf_z[ii][jj][kk] - ez_2_nl[ii][jj][kk - 1] * epsinf_z[ii][jj][kk - 1])) +
														 lat.inv_d[0] * (px_cp_2_nl[ii][jj][kk] - px_cp_2_nl[ii - 1][jj][kk]) +
														 lat.inv_d[1] * (py_cp_2_nl[ii][jj][kk] - py_cp_2_nl[ii][jj - 1][kk]) +
														 lat.inv_d[2] * (pz_cp_2_nl[ii][jj][kk] - pz_cp_2_nl[ii][jj][kk - 1]));

						nd_z_prev[ii][jj][kk] = (-1.0 / E0) * (EP0 * (lat.inv_d[0] * (ex_2_nl[ii][jj][kk] - ex_2_nl[ii - 1][jj][kk]) +
																	 lat.inv_d[1] * (ey_2_nl[ii][jj][kk] - ey_2_nl[ii][jj - 1][kk]) +
																	 lat.inv_d[2] * (ez_2_nl[ii][jj][kk] - ez_2_nl[ii][jj][kk - 1])));
					}
				}
			}
		}
	}

	fdtd::grid::exchange_margin(nd_x);
	fdtd::grid::exchange_margin(nd_z);

	if (fdtd::materials::model_spec().harmonic)
	{
		fdtd::grid::exchange_margin(nd_y);
		fdtd::grid::exchange_margin(nd_z_prev);
	}
}

/**
 * @brief Function that calls hudrodynamic update functions
 *
 * @details The hydrodynamic fields are updated using a 2-step process
 * see documentation. fdtd::materials::hydro::update_px() gives us an approximation for jx^{n+1/2}
 * and fdtd::materials::hydro::update_px_2() adds the nonlinear correction. A lot of fdtd::grid::exchange_margins()
 * functions are called to ensure each process has the nessecary field data at
 * it's boundaries.
 *
 */
void fdtd::materials::hydro::update_hydrodynamics()
{
	auto &ds = fdtd::materials::dispersive_state();
	int i, j, k;
	// Drude Update
	fdtd::materials::hydro::update_px();
	fdtd::materials::hydro::update_py();
	fdtd::materials::hydro::update_pz();

	// Switch here because we need all of the data from the process
	if (fdtd::materials::model_spec().antenna_material - DispMateOffsetD >= ModLorentz)
		fdtd::grid::exchange_margins(ds.c[0].j[fdtd::materials::model_spec().antenna_material - DispMateOffsetD], ds.c[1].j[fdtd::materials::model_spec().antenna_material - DispMateOffsetD], ds.c[2].j[fdtd::materials::model_spec().antenna_material - DispMateOffsetD]);
	if (fdtd::materials::model_spec().two_temperature)
		fdtd::grid::exchange_margins(ds.c[0].j[350 - DispMateOffsetD], ds.c[1].j[350 - DispMateOffsetD], ds.c[2].j[350 - DispMateOffsetD]);

	fdtd::materials::hydro::update_px_2();
	fdtd::materials::hydro::update_py_2();
	fdtd::materials::hydro::update_pz_2();

	// No need to switch here

  for (i = 0; i <= subdomain.cells_x + 2; i++) 
  {
    for (j = 0; j <= subdomain.cells_y + 2; j++) 
	{
      for (k = 0; k <= subdomain.cells_z + 2; k++) 
	  {
        if (i <= subdomain.cells_x && j <= subdomain.cells_y && k <= subdomain.cells_z) 
		{
          hx_prev[i][j][k] = Hx[i][j][k];
          hy_prev[i][j][k] = Hy[i][j][k];
          hz_prev[i][j][k] = Hz[i][j][k];
        }
        if (!fdtd::materials::model_spec().harmonic) 
		{
          if (fdtd::materials::model_spec().two_temperature) 
		  {
            ds.c[0].j[350 - DispMateOffsetD][i][j][k] += jx_extra[350 - DispMateOffsetD][i][j][k];
            ds.c[1].j[350 - DispMateOffsetD][i][j][k] += jy_extra[350 - DispMateOffsetD][i][j][k];
            ds.c[2].j[350 - DispMateOffsetD][i][j][k] += jz_extra[350 - DispMateOffsetD][i][j][k];
          }

					if (fdtd::materials::model_spec().antenna_material - DispMateOffsetD >= ModLorentz)
					{
						ds.c[0].j[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k] += jx_extra[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k];
						ds.c[1].j[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k] += jy_extra[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k];
						ds.c[2].j[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k] += jz_extra[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k];
					}
				}
				else
				{ // IF perturbative and therefore no TTM
					if (fdtd::materials::model_spec().antenna_material - DispMateOffsetD >= ModLorentz)
					{
						ds.c[0].j_nl[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k] += jx_extra[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k];
						ds.c[1].j_nl[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k] += jy_extra[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k];
						ds.c[2].j_nl[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k] += jz_extra[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k];
					}
				}
				if (fdtd::materials::model_spec().antenna_material - DispMateOffsetD >= ModLorentz)
				{
					jx_prev[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k] = ds.c[0].j[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k];
					jy_prev[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k] = ds.c[1].j[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k];
					jz_prev[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k] = ds.c[2].j[fdtd::materials::model_spec().antenna_material - DispMateOffsetD][i][j][k];
				}
				if (fdtd::materials::model_spec().two_temperature)
				{
					jx_prev[350 - DispMateOffsetD][i][j][k] = ds.c[0].j[350 - DispMateOffsetD][i][j][k];
					jy_prev[350 - DispMateOffsetD][i][j][k] = ds.c[1].j[350 - DispMateOffsetD][i][j][k];
					jz_prev[350 - DispMateOffsetD][i][j][k] = ds.c[2].j[350 - DispMateOffsetD][i][j][k];
				}
			}
		}
	}

	if (fdtd::materials::model_spec().harmonic)
	{
		fdtd::grid::exchange_margins(ds.c[0].j_nl[fdtd::materials::model_spec().antenna_material - DispMateOffsetD], ds.c[1].j_nl[fdtd::materials::model_spec().antenna_material - DispMateOffsetD], ds.c[2].j_nl[fdtd::materials::model_spec().antenna_material - DispMateOffsetD]);
	}
	if (fdtd::materials::model_spec().antenna_material - DispMateOffsetD >= ModLorentz)
		fdtd::grid::exchange_margins(jx_prev[fdtd::materials::model_spec().antenna_material - DispMateOffsetD], jy_prev[fdtd::materials::model_spec().antenna_material - DispMateOffsetD], jz_prev[fdtd::materials::model_spec().antenna_material - DispMateOffsetD]);
	if (fdtd::materials::model_spec().two_temperature)
		fdtd::grid::exchange_margins(jx_prev[350 - DispMateOffsetD], jy_prev[350 - DispMateOffsetD], jz_prev[350 - DispMateOffsetD]);
}

/**
 * @brief Function that populates epsinf_x
 *
 * @details Uses the same structure as the update functions to find which cells
 * countain a hydrodyanmic material, and then populates epsinf_x with epre_inf for
 * the current hydrodynamic material. Elsewhere, it is set to epre. This
 * function needs to be called before any simulation require the electron density
 *
 */

void fdtd::materials::hydro::setup_eps_inf_x() {
  const auto &mp = fdtd::materials::model_parameters();
  int i, j, k, u, v, w, U, V, W, mx, mmx; //
  int ii, jj, kk;
  for (i = 0; i < subdomain.cells_x; ++i) {
    ii = i + 1;

		if (i >= subdomain.i_cpml_min && i < subdomain.i_cpml_max)
		{

			U = 1;
			u = i - subdomain.inner_org_x;
		}
		else
		{

			if (i < subdomain.i_cpml_min)
			{

				U = 0;
			}
			else
			{

				U = 2;
			}
		}

		for (j = 1; j <= subdomain.owned_cells_y; ++j)
		{
			jj = j + 1;

			if (j >= subdomain.j_cpml_min && j <= subdomain.j_cpml_max)
			{

				V = 1;
				v = j - subdomain.inner_org_y;
			}
			else
			{

				if (j < subdomain.j_cpml_min)
				{

					V = 0;
				}
				else
				{

					V = 2;
				}

			}

			for (k = 1; k <= subdomain.owned_cells_z; ++k)
			{
				kk = k + 1;

				epsinf_x[ii][jj][kk] = 1.0;
				if (k >= subdomain.k_cpml_min && k <= subdomain.k_cpml_max)
				{
					W = 1;
					w = k - subdomain.inner_org_z;
				}
				else
				{

					if (k < subdomain.k_cpml_min)
					{

						W = 0;
					}
					else
					{

						W = 2;
					}

				}

				if ((U == 1) && (V == 1) && (W == 1))
				{
					mx = int(matind_x[u][v][w]);
				}
				else
				{

					if (V == 2)
					{
						mx = fdtd::materials::model_spec().medium_second;
					}
					else if (V == 0)
					{
						mx = fdtd::materials::model_spec().medium_first;
					}
					else
					{
						mx = background_y[j + subdomain.j_for_uni];
					}
				}
				if (mx >= DispMateOffsetD && mx < NRMAXMATERIALS)
				{
					mmx = mx - DispMateOffsetD;
					if (mmx >= ModLorentz && mmx < NonLinFreeElectron)
					{

						epsinf_x[ii][jj][kk] = mp.epre_inf[mmx];
					}

					if (mmx >= NonLinFreeElectron && mmx < NonlinearDrudeLorentz)
					{

						epsinf_x[ii][jj][kk] = mp.epre_inf[mmx];
					}
				}
				else
				{
					epsinf_x[ii][jj][kk] = mp.epre[mx];
				}
			}
		}
	}

	return;
}

/**
 * @brief Function that populates epsinf_y
 *
 * @details Uses the same structure as the update functions to find which cells
 * countain a hydrodyanmic material, and then populates epsinf_y with epre_inf for
 * the current hydrodynamic material. Elsewhere, it is set to epre. This
 * function needs to be called before any simulation require the electron density
 *
 */
void fdtd::materials::hydro::setup_eps_inf_y()
{
	const auto &mp = fdtd::materials::model_parameters();
	int i, j, k, u, v, w, U, V, W, my, mmy; //
	int ii, jj, kk;

	for (i = 1; i <= subdomain.owned_cells_x; ++i)
	{
		ii = i + 1;

		if (i >= subdomain.i_cpml_min && i <= subdomain.i_cpml_max)
		{

			U = 1;
			u = i - subdomain.inner_org_x;
		}
		else
		{

			if(i < subdomain.i_cpml_min) 
			{
				U = 0;
			}
			else 
			{
			U = 2;
			}
		}

		for (j = 0; j < subdomain.cells_y; ++j) {
			jj = j + 1;
			if (j >= subdomain.j_cpml_min && j < subdomain.j_cpml_max) {
				V = 1;
				v = j - subdomain.inner_org_y;
			}
			else {
				if (j < subdomain.j_cpml_min) {
					V = 0;
				}
				else {
					V = 2;
				}
			}

			for (k = 1; k <= subdomain.owned_cells_z; ++k)
			{
				kk = k + 1;
				epsinf_y[ii][jj][kk] = 1.0;

				if (k >= subdomain.k_cpml_min && k <= subdomain.k_cpml_max)
				{

					W = 1;
					w = k - subdomain.inner_org_z;
				}
				else
				{

					if (k < subdomain.k_cpml_min)
					{

						W = 0;
					}
					else
					{

						W = 2;
					}

				}

				if ((U == 1) && (V == 1) && (W == 1))
				{
					my = int(matind_y[u][v][w]);
				}
				else
				{

					if (V == 2)
						my = fdtd::materials::model_spec().medium_second;

					else if (V == 0)
						my = fdtd::materials::model_spec().medium_first;

					else
						my = background_y[j + subdomain.j_for_uni];
				}
				if (my >= DispMateOffsetD && my < NRMAXMATERIALS)
				{
					mmy = my - DispMateOffsetD;
					if (mmy >= ModLorentz && mmy < NonLinFreeElectron)
					{ // Free Electron Hydrodynamics model

						epsinf_y[ii][jj][kk] = mp.epre_inf[mmy];
					}

					if (mmy >= NonLinFreeElectron && mmy < NonlinearDrudeLorentz)
					{

						epsinf_y[ii][jj][kk] = mp.epre_inf[mmy];
					}
				}
				else
				{
					epsinf_y[ii][jj][kk] = mp.epre[my];
				}
			}
		}
	}

	return;
}

/**
 * @brief Function that populates epsinf_z
 *
 * @details Uses the same structure as the update functions to find which cells
 * countain a hydrodyanmic material, and then populates epsinf_z with epre_inf for
 * the current hydrodynamic material. Elsewhere, it is set to epre. This
 * function needs to be called before any simulation require the electron density
 *
 */
void fdtd::materials::hydro::setup_eps_inf_z()
{
	const auto &mp = fdtd::materials::model_parameters();
	int i, j, k, u, v, w, U, V, W, mz, mmz; //
	int ii, jj, kk;

	for (i = 1; i <= subdomain.owned_cells_x; ++i)
	{
		ii = i + 1;

		if (i >= subdomain.i_cpml_min && i <= subdomain.i_cpml_max)
		{

			U = 1;
			u = i - subdomain.inner_org_x;
		}
		else
		{

			if (i < subdomain.i_cpml_min)
			{

				U = 0;
			}
			else
			{

				U = 2;
			}

		}

		for (j = 1; j <= subdomain.owned_cells_y; ++j)
		{
			jj = j + 1;

			if (j >= subdomain.j_cpml_min && j <= subdomain.j_cpml_max)
			{
				V = 1;
				v = j - subdomain.inner_org_y;
			}
			else
			{

				if (j < subdomain.j_cpml_min)
				{

					V = 0;
				}
				else
				{

					V = 2;
				}

			}

			for (k = 0; k < subdomain.cells_z; ++k) {
				kk = k + 1;
				epsinf_z[ii][jj][kk] = 1.0;

				if (k >= subdomain.k_cpml_min && k < subdomain.k_cpml_max)
				{

					W = 1;
					w = k - subdomain.inner_org_z;
				}
				else
				{

					if (k < subdomain.k_cpml_min)
					{
         
						W = 0;
					}
					else
					{

						W = 2;
					}
				}

				if ((U == 1) && (V == 1) && (W == 1))
				{
					mz = int(matind_z[u][v][w]);
				}
				else
				{

					if (V == 2)
						mz = fdtd::materials::model_spec().medium_second;

					else if (V == 0)
						mz = fdtd::materials::model_spec().medium_first;

					else
						mz = background_y[j + subdomain.j_for_uni];
				}
				if (mz >= DispMateOffsetD && mz < NRMAXMATERIALS)
				{
					mmz = mz - DispMateOffsetD;

					if (mmz >= ModLorentz && mmz < NonLinFreeElectron)
					{ // Free Electron Hydrodynamics model
						epsinf_z[ii][jj][kk] = mp.epre_inf[mmz];
					}

					if (mmz >= NonLinFreeElectron && mmz < NonlinearDrudeLorentz)
					{

						epsinf_z[ii][jj][kk] = mp.epre_inf[mmz];
					}
				}
				else
				{
					epsinf_z[ii][jj][kk] = mp.epre[mz];
				}
			}
		}
	}
	return;
}

/**
 * @brief 1st step of update algorithm for hydrodynamic model + calculation of
 * MULTIPOLES, only valid for small particles
 *
 * @details Updates jx using a "linear/local" approximation - the first step in
 * the hydrodyanmic update procedure. Also supports the nonlinear Drude-Lorentz
 * model + TTM. Uniquely, this function also calculates the
 * electron density integrals required for calculating the MULTIPOLES of a
 * nanoparticle, who's diameter << wavelength.
 *
 * @note The multipole calculation makes this function more confusing to look at.
 * Probably best to remove this.
 *
 */

void fdtd::materials::hydro::update_px() {
  auto &ds = fdtd::materials::dispersive_state();
  int i, j, k, u, v, w, U, V, W, mx, mmx; //
  int ii, jj, kk;
  double Grad_Div, Grad_Div2;
  double ND1, powerapprox;
  double gamma;

  for (i = 0; i < subdomain.cells_x; ++i) {
    ii = i + 1;

    if (i >= subdomain.i_cpml_min && i < subdomain.i_cpml_max) {

      U = 1;
      u = i - subdomain.inner_org_x;

    } else {

      if (i < subdomain.i_cpml_min) {

        U = 0;

      } else {

        U = 2;
      }
    }

    for (j = 1; j <= subdomain.owned_cells_y; ++j) {
      jj = j + 1;

      if (j >= subdomain.j_cpml_min && j <= subdomain.j_cpml_max) {

        V = 1;
        v = j - subdomain.inner_org_y;

      } else {

        if (j < subdomain.j_cpml_min) {

          V = 0;

        } else {

          V = 2;
        }

      }

      for (k = 1; k <= subdomain.owned_cells_z; ++k) {
        kk = k + 1;

        if (k >= subdomain.k_cpml_min && k <= subdomain.k_cpml_max) {
          W = 1;
          w = k - subdomain.inner_org_z;

        } else {

          if (k < subdomain.k_cpml_min) {

            W = 0;

          } else {

            W = 2;
          }

        }

        if ((U == 1) && (V == 1) && (W == 1)) {
          mx = int(matind_x[u][v][w]);
        } else {

          if (V == 2) {
            mx = fdtd::materials::model_spec().medium_second;
          } else if (V == 0) {
            mx = fdtd::materials::model_spec().medium_first;
          } else {
            mx = background_y[j + subdomain.j_for_uni];

          }

        }

        if (mx >= DispMateOffsetD && mx < NRMAXMATERIALS) {
          mmx = mx - DispMateOffsetD;
          if (mmx >= ModLorentz && mmx < NonLinFreeElectron) {

            if (with_pressure) Grad_Div = lat.inv_d[0] * (nd_x[ii + 1][jj][kk] - nd_x[ii][jj][kk]);
            else Grad_Div = 0.0;
            ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii + 1][jj][kk]);

            if (fdtd::materials::model_spec().harmonic) {


              if (with_pressure) Grad_Div2 = lat.inv_d[0] * (nd_y[ii + 1][jj][kk] - nd_y[ii][jj][kk]);
              else Grad_Div2 = 0.0;
              powerapprox = powerterm[mmx] * (1.0 + (2.0 / 3.0) * ND1 * inverse_edens[mmx]);
              
              ds.c[0].j[mmx][ii][jj][kk] = fdtd::materials::electron_gas(mmx).nl1 * jx_prev[mmx][ii][jj][kk] + (fdtd::materials::electron_gas(mmx).nl2) * Ex[i][j][k] + fdtd::materials::electron_gas(mmx).nl5 * powerterm[mmx] * Grad_Div;
              ds.c[0].j_nl[mmx][ii][jj][kk] = fdtd::materials::electron_gas(mmx).nl1 * ds.c[0].j_nl[mmx][ii][jj][kk] + (fdtd::materials::electron_gas(mmx).nl2 * ex_nl[i][j][k] + fdtd::materials::electron_gas(mmx).nl3 * ND1 * Ex[i][j][k]) + fdtd::materials::electron_gas(mmx).nl5 * powerterm[mmx] * Grad_Div2 + fdtd::materials::electron_gas(mmx).nl5 * powerapprox * Grad_Div;

            } else {
              ds.c[0].j[mmx][ii][jj][kk] = fdtd::materials::electron_gas(mmx).nl1 * jx_prev[mmx][ii][jj][kk] + (fdtd::materials::electron_gas(mmx).nl2 + fdtd::materials::electron_gas(mmx).nl3 * ND1) * Ex[i][j][k];
            }

          }

          else if (mmx >= NonLinFreeElectron && mmx < NonlinearDrudeLorentz) {

            gamma = fdtd::materials::thermal::damping_constant(fdtd::materials::thermal::electron_temp[i][j][k], fdtd::materials::thermal::lattice_temp[i][j][k]);
            ds.c[0].j[mmx][ii][jj][kk] = ((1.0 / fdtd::materials::model_spec().dt - 0.5 * gamma) / (1.0 / fdtd::materials::model_spec().dt + 0.5 * gamma)) * jx_prev[mmx][ii][jj][kk] + (EP0 * fdtd::materials::thermal::plasma_frequency((fdtd::materials::thermal::electron_temp[i][j][k])) * fdtd::materials::thermal::plasma_frequency((fdtd::materials::thermal::electron_temp[i][j][k])) / (1.0 / fdtd::materials::model_spec().dt + 0.5 * gamma)) * Ex[i][j][k];
          }

        }
      }
    }
  }

  return;

}

/**
 * @brief 2nd step of update algorithm for hydrodynamic model
 *
 * @details Calculates the nonlinear correction to jx, see docs
 */

void fdtd::materials::hydro::update_px_2() {
  auto &ds = fdtd::materials::dispersive_state();
  int i, j, k, u, v, w, U, V, W, mx, mmx; //
  int ii, jj, kk;
  double Grad_Div;
  double GradNY, GradNZ;
  double Vx1, Vy1, Vz1, Hz1, Hy1, Hz2, Hy2, VdotGrad, VcrossH, ND1, DivJ, VdotGradN, VdotGradNeq;
  double Px_this, Py_this, Pz_this, Px_this2, Py_this2, Pz_this2, Px_thisxp1, Py_thisxp1, Pz_thisxp1, Px_this2xp1, Py_this2xp1, Pz_this2xp1;
  double Px_thisxm1, Px_this2xm1, Py_thisym1, Py_this2ym1, Py_thisym1xp1, Py_this2ym1xp1, Pz_thiszm1, Pz_this2zm1, Pz_thiszm1xp1, Pz_this2zm1xp1, TOTND1;
  double NONLINDrude4TTM;
  for (i = 0; i < subdomain.cells_x; ++i) {
    ii = i + 1;

    if (i >= subdomain.i_cpml_min && i < subdomain.i_cpml_max) {

      U = 1;
      u = i - subdomain.inner_org_x;

    } else {

      if (i < subdomain.i_cpml_min) {

        U = 0;

      } else {

        U = 2;
      }
    }

    for (j = 1; j <= subdomain.owned_cells_y; ++j) {
      jj = j + 1;

      if (j >= subdomain.j_cpml_min && j <= subdomain.j_cpml_max) {

        V = 1;
        v = j - subdomain.inner_org_y;

      } else {

        if (j < subdomain.j_cpml_min) {

          V = 0;

        } else {

          V = 2;
        }

      }

      for (k = 1; k <= subdomain.owned_cells_z; ++k) {
        kk = k + 1;

        if (k >= subdomain.k_cpml_min && k <= subdomain.k_cpml_max) {

          W = 1;
          w = k - subdomain.inner_org_z;

        } else {

          if (k < subdomain.k_cpml_min) {

            W = 0;

          } else {

            W = 2;
          }

        }

        if ((U == 1) && (V == 1) && (W == 1)) {
          mx = int(matind_x[u][v][w]);

        } else {

          if (V == 2) {
            mx = fdtd::materials::model_spec().medium_second;
          } else if (V == 0) {
            mx = fdtd::materials::model_spec().medium_first;
          } else {
            mx = background_y[j + subdomain.j_for_uni];

          }

        }
        if (mx >= DispMateOffsetD && mx < NRMAXMATERIALS) {
          mmx = mx - DispMateOffsetD;
          if (mmx >= ModLorentz && mmx < NonLinFreeElectron) {

            if (use_backwards_differences) {
              if (with_pressure) Grad_Div = lat.inv_d[0] * (nd_x[ii + 1][jj][kk] - nd_x[ii][jj][kk]);
              else Grad_Div = 0.0;
              ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii + 1][jj][kk]);
              TOTND1 = ND1 + equilibrium_electron_density[mmx];

              Px_this = jx_prev[mmx][ii][jj][kk];
              Py_this = jy_prev[mmx][ii][jj][kk];
              Pz_this = jz_prev[mmx][ii][jj][kk];

              Px_thisxp1 = jx_prev[mmx][ii + 1][jj][kk];
              Py_thisxp1 = jy_prev[mmx][ii + 1][jj][kk];
              Pz_thisxp1 = jz_prev[mmx][ii + 1][jj][kk];

              Px_thisxm1 = jx_prev[mmx][ii - 1][jj][kk];

              Py_thisym1 = jy_prev[mmx][ii][jj - 1][kk];
              Py_thisym1xp1 = jy_prev[mmx][ii + 1][jj - 1][kk];

              Pz_thiszm1 = jz_prev[mmx][ii][jj][kk - 1];
              Pz_thiszm1xp1 = jz_prev[mmx][ii + 1][jj][kk - 1];

              Vx1 = (Px_this);
              Vy1 = (Py_this);
              Vz1 = (Pz_this);

              Hy2 = (hy_prev[i][j][k]);
              Hz2 = (hz_prev[i][j][k]);
              Hy1 = (Hy[i][j][k]);
              Hz1 = (Hz[i][j][k]);

              Hy1 = 0.5 * (Hy1 + Hy2);
              Hz1 = 0.5 * (Hz1 + Hz2);

              if (with_convection) {
                VdotGrad = (Vx1 * (Px_this - Px_thisxm1) * lat.inv_d[0] +
                  Vy1 * (Px_this - jx_prev[mmx][ii][jj - 1][kk]) * lat.inv_d[1] +
                  Vz1 * (Px_this - jx_prev[mmx][ii][jj][kk - 1]) * lat.inv_d[2]);

              } else VdotGrad = 0.0;

              if (with_mag_field) {
                VcrossH = Vy1 * Hz1 - Vz1 * Hy1;
              } else VcrossH = 0.0;

              DivJ = (lat.inv_d[0] * (Px_this - Px_thisxm1) +
                lat.inv_d[1] * (Py_this - Py_thisym1) +
                lat.inv_d[2] * (Pz_this - Pz_thiszm1));

              if (!fdtd::materials::model_spec().harmonic) {
                GradNY = (nd_x[ii][jj][kk] - nd_x[ii][jj - 1][kk]) * lat.inv_d[1];
                GradNZ = (nd_x[ii][jj][kk] - nd_x[ii][jj][kk - 1]) * lat.inv_d[2];

                VdotGradN = Vx1 * Grad_Div + Vy1 * GradNY + Vz1 * GradNZ;
              }
              VdotGradN = 0.0;
              
              VdotGradNeq = 0.0;

              if (!fdtd::materials::model_spec().harmonic) {
                jx_extra[mmx][ii][jj][kk] = fdtd::materials::electron_gas(mmx).nl4 * (1.0 / (E0 * (TOTND1))) * (VdotGrad + Vx1 * (DivJ -
                    (VdotGradN + VdotGradNeq) / (TOTND1))) - fdtd::materials::electron_gas(mmx).nl4 * (MU0 * E0 / effective_electron_mass[mmx]) * VcrossH +
                  fdtd::materials::electron_gas(mmx).nl5 * pow(TOTND1, 2.0 / 3.0) * Grad_Div;
              } else {
                jx_extra[mmx][ii][jj][kk] = fdtd::materials::electron_gas(mmx).nl4 * inverse_chargedens[mmx] * (VdotGrad + Vx1 * (DivJ)) -
                  fdtd::materials::electron_gas(mmx).nl4 * (MU0 * E0 / effective_electron_mass[mmx]) * VcrossH;
              }
            } else {
              if (with_pressure) Grad_Div = lat.inv_d[0] * (nd_x[ii + 1][jj][kk] - nd_x[ii][jj][kk]);
              else Grad_Div = 0.0;
              ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii + 1][jj][kk]);
              TOTND1 = ND1 + equilibrium_electron_density[mmx];

              Px_this = jx_prev[mmx][ii][jj][kk];
              Py_this = jy_prev[mmx][ii][jj][kk];
              Pz_this = jz_prev[mmx][ii][jj][kk];
              Px_this2 = ds.c[0].j[mmx][ii][jj][kk];
              Py_this2 = ds.c[1].j[mmx][ii][jj][kk];
              Pz_this2 = ds.c[2].j[mmx][ii][jj][kk];
              Px_thisxp1 = jx_prev[mmx][ii + 1][jj][kk];
              Py_thisxp1 = jy_prev[mmx][ii + 1][jj][kk];
              Pz_thisxp1 = jz_prev[mmx][ii + 1][jj][kk];
              Px_this2xp1 = ds.c[0].j[mmx][ii + 1][jj][kk];
              Py_this2xp1 = ds.c[1].j[mmx][ii + 1][jj][kk];
              Pz_this2xp1 = ds.c[2].j[mmx][ii + 1][jj][kk];

              Px_thisxm1 = jx_prev[mmx][ii - 1][jj][kk];
              Px_this2xm1 = ds.c[0].j[mmx][ii - 1][jj][kk];

              Py_thisym1 = jy_prev[mmx][ii][jj - 1][kk];
              Py_this2ym1 = ds.c[1].j[mmx][ii][jj - 1][kk];
              Py_thisym1xp1 = jy_prev[mmx][ii + 1][jj - 1][kk];
              Py_this2ym1xp1 = ds.c[1].j[mmx][ii + 1][jj - 1][kk];

              Pz_thiszm1 = jz_prev[mmx][ii][jj][kk - 1];
              Pz_this2zm1 = ds.c[2].j[mmx][ii][jj][kk - 1];
              Pz_thiszm1xp1 = jz_prev[mmx][ii + 1][jj][kk - 1];
              Pz_this2zm1xp1 = ds.c[2].j[mmx][ii + 1][jj][kk - 1];

              Vx1 = (Px_this + Px_this2) * 0.5;

              Vy1 = 0.25 * (Py_this + Py_thisxp1 + Py_thisym1 + Py_thisym1xp1);
              Vy1 = Vy1 + 0.25 * (Py_this2 + Py_this2xp1 + Py_this2ym1 + Py_this2ym1xp1);
              Vy1 /= 2.0;

              Vz1 = 0.25 * (Pz_this + Pz_thisxp1 + Pz_thiszm1 + Pz_thiszm1xp1);
              Vz1 = Vz1 + 0.25 * (Pz_this2 + Pz_this2xp1 + Pz_this2zm1 + Pz_this2zm1xp1);
              Vz1 /= 2.0;

              Hy2 = 0.5 * (hy_prev[i][j][k] + hy_prev[i][j][k - 1]);
              Hz2 = 0.5 * (hz_prev[i][j][k] + hz_prev[i][j - 1][k]);
              Hy1 = 0.5 * (Hy[i][j][k] + Hy[i][j][k - 1]);
              Hz1 = 0.5 * (Hz[i][j][k] + Hz[i][j - 1][k]);

              Hy1 = 0.5 * (Hy1 + Hy2);
              Hz1 = 0.5 * (Hz1 + Hz2);

              if (with_convection) {

                VdotGrad = 0.25 * (Vx1 * (Px_thisxp1 + Px_this2xp1 - Px_thisxm1 - Px_this2xm1) * lat.inv_d[0] +
                  Vy1 * (jx_prev[mmx][ii][jj + 1][kk] + ds.c[0].j[mmx][ii][jj + 1][kk] - jx_prev[mmx][ii][jj - 1][kk] - ds.c[0].j[mmx][ii][jj - 1][kk]) * lat.inv_d[1] +
                  Vz1 * (jx_prev[mmx][ii][jj][kk + 1] + ds.c[0].j[mmx][ii][jj][kk + 1] - jx_prev[mmx][ii][jj][kk - 1] - ds.c[0].j[mmx][ii][jj][kk - 1]) * lat.inv_d[2]);

              } else VdotGrad = 0.0;

              if (with_mag_field) {
                VcrossH = Vy1 * Hz1 - Vz1 * Hy1;
              } else VcrossH = 0.0;

              DivJ = 0.25 * (lat.inv_d[0] * (Px_thisxp1 + Px_this2xp1 - Px_thisxm1 - Px_this2xm1) +
                lat.inv_d[1] * (Py_this + Py_thisxp1 + Py_this2 + Py_this2xp1 -
                  Py_thisym1 - Py_thisym1xp1 - Py_this2ym1 - Py_this2ym1xp1) +
                lat.inv_d[2] * (Pz_this + Pz_thisxp1 + Pz_this2 + Pz_this2xp1 -
                  Pz_thiszm1 - Pz_thiszm1xp1 - Pz_this2zm1 - Pz_this2zm1xp1));

              if (!fdtd::materials::model_spec().harmonic) {
                GradNY = (nd_x[ii + 1][jj + 1][kk] + nd_x[ii][jj + 1][kk] - nd_x[ii + 1][jj - 1][kk] - nd_x[ii][jj - 1][kk]) * 0.25 * lat.inv_d[1];
                GradNZ = (nd_x[ii + 1][jj][kk + 1] + nd_x[ii][jj][kk + 1] - nd_x[ii + 1][jj][kk - 1] - nd_x[ii][jj][kk - 1]) * 0.25 * lat.inv_d[2];

                VdotGradN = Vx1 * Grad_Div + Vy1 * GradNY + Vz1 * GradNZ;
              }

              VdotGradNeq = 0.0;

              if (!fdtd::materials::model_spec().harmonic) {
                jx_extra[mmx][ii][jj][kk] = fdtd::materials::electron_gas(mmx).nl4 * (1.0 / (E0 * (TOTND1))) * (VdotGrad + Vx1 * (DivJ -
                  (VdotGradN + VdotGradNeq) / (TOTND1))) - fdtd::materials::electron_gas(mmx).nl4 * (MU0 * E0 / effective_electron_mass[mmx]) * VcrossH; // + fdtd::materials::electron_gas(mmx).nl5*pow(TOTND1,2.0/3.0)*Grad_Div;

              } else {
                jx_extra[mmx][ii][jj][kk] = fdtd::materials::electron_gas(mmx).nl4 * inverse_chargedens[mmx] * (VdotGrad + Vx1 * (DivJ)) -
                  fdtd::materials::electron_gas(mmx).nl4 * (MU0 * E0 / effective_electron_mass[mmx]) * VcrossH;
              }
            }
          }
        
          else if (mmx >= NonLinFreeElectron && mmx < NonlinearDrudeLorentz) {
            if (with_pressure) Grad_Div = lat.inv_d[0] * (nd_x[ii + 1][jj][kk] - nd_x[ii][jj][kk]);
            else Grad_Div = 0.0;
            ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii + 1][jj][kk]);
            TOTND1 = ND1 + equilibrium_electron_density[mmx];

            Px_this = jx_prev[mmx][ii][jj][kk];
            Py_this = jy_prev[mmx][ii][jj][kk];
            Pz_this = jz_prev[mmx][ii][jj][kk];
            Px_this2 = ds.c[0].j[mmx][ii][jj][kk];
            Py_this2 = ds.c[1].j[mmx][ii][jj][kk];
            Pz_this2 = ds.c[2].j[mmx][ii][jj][kk];
            Px_thisxp1 = jx_prev[mmx][ii + 1][jj][kk];
            Py_thisxp1 = jy_prev[mmx][ii + 1][jj][kk];
            Pz_thisxp1 = jz_prev[mmx][ii + 1][jj][kk];
            Px_this2xp1 = ds.c[0].j[mmx][ii + 1][jj][kk];
            Py_this2xp1 = ds.c[1].j[mmx][ii + 1][jj][kk];
            Pz_this2xp1 = ds.c[2].j[mmx][ii + 1][jj][kk];

            Px_thisxm1 = jx_prev[mmx][ii - 1][jj][kk];
            Px_this2xm1 = ds.c[0].j[mmx][ii - 1][jj][kk];

            Py_thisym1 = jy_prev[mmx][ii][jj - 1][kk];
            Py_this2ym1 = ds.c[1].j[mmx][ii][jj - 1][kk];
            Py_thisym1xp1 = jy_prev[mmx][ii + 1][jj - 1][kk];
            Py_this2ym1xp1 = ds.c[1].j[mmx][ii + 1][jj - 1][kk];

            Pz_thiszm1 = jz_prev[mmx][ii][jj][kk - 1];
            Pz_this2zm1 = ds.c[2].j[mmx][ii][jj][kk - 1];
            Pz_thiszm1xp1 = jz_prev[mmx][ii + 1][jj][kk - 1];
            Pz_this2zm1xp1 = ds.c[2].j[mmx][ii + 1][jj][kk - 1];

            Vx1 = (Px_this + Px_this2) * 0.5;

            Vy1 = 0.25 * (Py_this + Py_thisxp1 + Py_thisym1 + Py_thisym1xp1);
            Vy1 = Vy1 + 0.25 * (Py_this2 + Py_this2xp1 + Py_this2ym1 + Py_this2ym1xp1);
            Vy1 /= 2.0;

            Vz1 = 0.25 * (Pz_this + Pz_thisxp1 + Pz_thiszm1 + Pz_thiszm1xp1);
            Vz1 = Vz1 + 0.25 * (Pz_this2 + Pz_this2xp1 + Pz_this2zm1 + Pz_this2zm1xp1);
            Vz1 /= 2.0;

            Hy2 = 0.5 * (hy_prev[i][j][k] + hy_prev[i][j][k - 1]);
            Hz2 = 0.5 * (hz_prev[i][j][k] + hz_prev[i][j - 1][k]);
            Hy1 = 0.5 * (Hy[i][j][k] + Hy[i][j][k - 1]);
            Hz1 = 0.5 * (Hz[i][j][k] + Hz[i][j - 1][k]);

            Hy1 = 0.5 * (Hy1 + Hy2);
            Hz1 = 0.5 * (Hz1 + Hz2);

            if (with_convection) {

              VdotGrad = 0.25 * (Vx1 * (Px_thisxp1 + Px_this2xp1 - Px_thisxm1 - Px_this2xm1) * lat.inv_d[0] +
                Vy1 * (jx_prev[mmx][ii][jj + 1][kk] + ds.c[0].j[mmx][ii][jj + 1][kk] - jx_prev[mmx][ii][jj - 1][kk] - ds.c[0].j[mmx][ii][jj - 1][kk]) * lat.inv_d[1] +
                Vz1 * (jx_prev[mmx][ii][jj][kk + 1] + ds.c[0].j[mmx][ii][jj][kk + 1] - jx_prev[mmx][ii][jj][kk - 1] - ds.c[0].j[mmx][ii][jj][kk - 1]) * lat.inv_d[2]);

            } else VdotGrad = 0.0;

            if (with_mag_field) {
              VcrossH = Vy1 * Hz1 - Vz1 * Hy1;
            } else VcrossH = 0.0;

            DivJ = 0.25 * (lat.inv_d[0] * (Px_thisxp1 + Px_this2xp1 - Px_thisxm1 - Px_this2xm1) +
              lat.inv_d[1] * (Py_this + Py_thisxp1 + Py_this2 + Py_this2xp1 -
                Py_thisym1 - Py_thisym1xp1 - Py_this2ym1 - Py_this2ym1xp1) +
              lat.inv_d[2] * (Pz_this + Pz_thisxp1 + Pz_this2 + Pz_this2xp1 -
                Pz_thiszm1 - Pz_thiszm1xp1 - Pz_this2zm1 - Pz_this2zm1xp1));

            VdotGradN = 0.0;
            VdotGradNeq = 0.0;

            NONLINDrude4TTM = 1.0 / (1.0 / fdtd::materials::model_spec().dt + 0.5 * fdtd::materials::thermal::damping_constant(fdtd::materials::thermal::electron_temp[i][j][k], fdtd::materials::thermal::lattice_temp[i][j][k]));

            jx_extra[mmx][ii][jj][kk] = NONLINDrude4TTM * (1.0 / (E0 * (TOTND1))) * (VdotGrad + Vx1 * (DivJ -
              (VdotGradN + VdotGradNeq) / (TOTND1))) - NONLINDrude4TTM * ((MU0 * EP0 * fdtd::materials::thermal::plasma_frequency(fdtd::materials::thermal::electron_temp[i][j][k]) * fdtd::materials::thermal::plasma_frequency(fdtd::materials::thermal::electron_temp[i][j][k])) / (E0 * TOTND1)) * VcrossH; // + fdtd::materials::electron_gas(mmx).nl5*pow(TOTND1,2.0/3.0)*Grad_Div;

          }

        }
      }
    }
  }

  return;

}

/**
 * @brief 1st step of update algorithm for hydrodynamic model
 *
 * @details Updates jy using a "linear/local" approximation - the first step in
 * the hydrodyanmic update procedure. Also supports the nonlinear Drude-Lorentz
 * model + TTM.
 */
void fdtd::materials::hydro::update_py()
{
	auto &ds = fdtd::materials::dispersive_state();
	int i, j, k, u, v, w, U, V, W, my, mmy; //
	int ii, jj, kk;
	double Grad_Div, Grad_Div2;
	double ND1, powerapprox;

	double gamma;



	for (i = 1; i <= subdomain.owned_cells_x; ++i)
	{
		ii = i + 1;

		if (i >= subdomain.i_cpml_min && i <= subdomain.i_cpml_max)
		{

			U = 1;
			u = i - subdomain.inner_org_x;
		}
		else
		{

			if (i < subdomain.i_cpml_min)
			{

				U = 0;
			}
			else
			{

				U = 2;
			}

		}

		for (j = 0; j < subdomain.cells_y; ++j) {
			jj = j + 1;

			if (j >= subdomain.j_cpml_min && j < subdomain.j_cpml_max)
			{

				V = 1;
				v = j - subdomain.inner_org_y;
			}
			else
			{

				if (j < subdomain.j_cpml_min)
				{

					V = 0;
				}
				else
				{

					V = 2;
				}
			}

			for (k = 1; k <= subdomain.owned_cells_z; ++k)
			{
				kk = k + 1;

				if (k >= subdomain.k_cpml_min && k <= subdomain.k_cpml_max)
				{

					W = 1;
					w = k - subdomain.inner_org_z;
				}
				else
				{

					if (k < subdomain.k_cpml_min)
					{

						W = 0;
					}
					else
					{

						W = 2;
					}

				}

				if ((U == 1) && (V == 1) && (W == 1))
				{
					my = int(matind_y[u][v][w]);
				}
				else
				{

					if (V == 2)
						my = fdtd::materials::model_spec().medium_second;

					else if (V == 0)
						my = fdtd::materials::model_spec().medium_first;

					else
						my = background_y[j + subdomain.j_for_uni];
				}

				// if(my==nl_code){my=NL_mat; my_=nl_code;} //(SHG || THG) &&
				if (my >= DispMateOffsetD && my < NRMAXMATERIALS)
				{
					mmy = my - DispMateOffsetD;
					if (mmy >= ModLorentz && mmy < NonLinFreeElectron)
					{ // Free Electron Hydrodynamics model
						if (with_pressure)
							Grad_Div = lat.inv_d[1] * (nd_x[ii][jj + 1][kk] - nd_x[ii][jj][kk]);
						else
							Grad_Div = 0.0;
						ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii][jj + 1][kk]);
						if (!fdtd::materials::model_spec().harmonic)
						{
							ds.c[1].j[mmy][ii][jj][kk] = fdtd::materials::electron_gas(mmy).nl1 * jy_prev[mmy][ii][jj][kk] + (fdtd::materials::electron_gas(mmy).nl2 + fdtd::materials::electron_gas(mmy).nl3 * ND1) * Ey[i][j][k];
						}
						else
						{
							powerapprox = powerterm[mmy] * (1.0 + (2.0 / 3.0) * ND1 * inverse_edens[mmy]);

							if (with_pressure)
								Grad_Div2 = lat.inv_d[1] * (nd_y[ii][jj + 1][kk] - nd_y[ii][jj][kk]);
							else
								Grad_Div2 = 0.0;

							ds.c[1].j[mmy][ii][jj][kk] = fdtd::materials::electron_gas(mmy).nl1 * jy_prev[mmy][ii][jj][kk] + (fdtd::materials::electron_gas(mmy).nl2) * Ey[i][j][k] + fdtd::materials::electron_gas(mmy).nl5 * powerterm[mmy] * Grad_Div;
							ds.c[1].j_nl[mmy][ii][jj][kk] = fdtd::materials::electron_gas(mmy).nl1 * ds.c[1].j_nl[mmy][ii][jj][kk] + (fdtd::materials::electron_gas(mmy).nl2 * ey_nl[i][j][k] + fdtd::materials::electron_gas(mmy).nl3 * ND1 * Ey[i][j][k]) + fdtd::materials::electron_gas(mmy).nl5 * powerterm[mmy] * Grad_Div2 + fdtd::materials::electron_gas(mmy).nl5 * powerapprox * Grad_Div;
						}
					}
					else if (mmy >= NonLinFreeElectron && mmy < NonlinearDrudeLorentz)
					{

						gamma = fdtd::materials::thermal::damping_constant(fdtd::materials::thermal::electron_temp[i][j][k], fdtd::materials::thermal::lattice_temp[i][j][k]);
						ds.c[1].j[mmy][ii][jj][kk] = ((1.0 / fdtd::materials::model_spec().dt - 0.5 * gamma) / (1.0 / fdtd::materials::model_spec().dt + 0.5 * gamma)) * jy_prev[mmy][ii][jj][kk] + (EP0 * fdtd::materials::thermal::plasma_frequency((fdtd::materials::thermal::electron_temp[i][j][k])) * fdtd::materials::thermal::plasma_frequency((fdtd::materials::thermal::electron_temp[i][j][k])) / (1.0 / fdtd::materials::model_spec().dt + 0.5 * gamma)) * Ey[i][j][k];
					}
				}
			}
		}
	}

	return;
}

/**
 * @brief 2nd step of update algorithm for hydrodynamic model
 *
 * @details Calculates the nonlinear correction to jy, see docs
 */
void fdtd::materials::hydro::update_py_2() {
  auto &ds = fdtd::materials::dispersive_state();
  int i, j, k, u, v, w, U, V, W, my, mmy; //
  int ii, jj, kk;
  double Grad_Div;
  double Vx1, Vy1, Vz1, Hx1, Hz1, Hx2, Hz2, VdotGrad, VcrossH, ND1, DivJ, VdotGradN, VdotGradNeq;
  double TOTND1;
  double GradNX, GradNZ;
  double Px_this, Py_this, Pz_this, Px_this2, Py_this2, Pz_this2, Px_thisyp1, Py_thisyp1, Pz_thisyp1, Px_this2yp1, Py_this2yp1, Pz_this2yp1, Py_thisym1, Py_this2ym1, Px_thisxm1, Px_this2xm1, Px_thisxm1yp1, Px_this2xm1yp1, Pz_thiszm1, Pz_this2zm1, Pz_thiszm1yp1, Pz_this2zm1yp1;
  double NONLINDrude4TTM;

  for (i = 1; i <= subdomain.owned_cells_x; ++i) {
    ii = i + 1;

    if (i >= subdomain.i_cpml_min && i <= subdomain.i_cpml_max) {

      U = 1;
      u = i - subdomain.inner_org_x;

    } else {

      if (i < subdomain.i_cpml_min) {

        U = 0;

      } else {

        U = 2;
      }

    }

    for (j = 0; j < subdomain.cells_y; ++j) {
      jj = j + 1;

      if (j >= subdomain.j_cpml_min && j < subdomain.j_cpml_max) {

        V = 1;
        v = j - subdomain.inner_org_y;

      } else {

        if (j < subdomain.j_cpml_min) {

          V = 0;

        } else {

          V = 2;
        }
      }

      for (k = 1; k <= subdomain.owned_cells_z; ++k) {
        kk = k + 1;

        if (k >= subdomain.k_cpml_min && k <= subdomain.k_cpml_max) {

          W = 1;
          w = k - subdomain.inner_org_z;

        } else {

          if (k < subdomain.k_cpml_min) {

            W = 0;

          } else {

            W = 2;
          }

        }

        if ((U == 1) && (V == 1) && (W == 1)) {
          my = int(matind_y[u][v][w]);
        } else {

          if (V == 2) my = fdtd::materials::model_spec().medium_second;

          else if (V == 0) my = fdtd::materials::model_spec().medium_first;

          else my = background_y[j + subdomain.j_for_uni];

        }
        
        if (my >= DispMateOffsetD && my < NRMAXMATERIALS) {
          mmy = my - DispMateOffsetD;
          if (mmy >= ModLorentz && mmy < NonLinFreeElectron) { //Free Electron Hydrodynamics model

            if (use_backwards_differences) {

              if (with_pressure) Grad_Div = lat.inv_d[1] * (nd_x[ii][jj + 1][kk] - nd_x[ii][jj][kk]);
              else Grad_Div = 0.0;
              ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii][jj + 1][kk]);
              TOTND1 = ND1 + equilibrium_electron_density[mmy];

              Px_this = jx_prev[mmy][ii][jj][kk];
              Py_this = jy_prev[mmy][ii][jj][kk];
              Pz_this = jz_prev[mmy][ii][jj][kk];
              Px_thisyp1 = jx_prev[mmy][ii][jj + 1][kk];
              Py_thisyp1 = jy_prev[mmy][ii][jj + 1][kk];
              Pz_thisyp1 = jz_prev[mmy][ii][jj + 1][kk];

              Py_thisym1 = jy_prev[mmy][ii][jj - 1][kk];

              Px_thisxm1 = jx_prev[mmy][ii - 1][jj][kk];
              Px_thisxm1yp1 = jx_prev[mmy][ii - 1][jj + 1][kk];

              Pz_thiszm1 = jz_prev[mmy][ii][jj][kk - 1];
              Pz_thiszm1yp1 = jz_prev[mmy][ii][jj + 1][kk - 1];

              Vy1 = (Py_this);
              Vx1 = (Px_this);
              Vz1 = (Pz_this);

              Hx1 = (Hx[i][j][k]);
              Hz1 = (Hz[i][j][k]);
              Hx2 = (hx_prev[i][j][k]);
              Hz2 = (hz_prev[i][j][k]);

              Hx1 = 0.5 * (Hx1 + Hx2);
              Hz1 = 0.5 * (Hz1 + Hz2);

              if (with_mag_field) {
                VcrossH = Vz1 * Hx1 - Vx1 * Hz1;
              } else VcrossH = 0.0;

              if (with_convection) {
                
                VdotGrad = (Vx1 * (Py_this - jy_prev[mmy][ii - 1][jj][kk]) * lat.inv_d[0] +
                  Vy1 * (Py_this - Py_thisym1) * lat.inv_d[1] +
                  Vz1 * (Py_this - jy_prev[mmy][ii][jj][kk - 1]) * lat.inv_d[2]);
              } else VdotGrad = 0.0;

              DivJ = (lat.inv_d[1] * (Py_this - Py_thisym1) +
                lat.inv_d[0] * (Px_this - Px_thisxm1) +
                lat.inv_d[2] * (Pz_this - Pz_thiszm1));

              VdotGradNeq = 0.0;

              if (!fdtd::materials::model_spec().harmonic) {
                GradNX = (nd_x[ii][jj][kk] - nd_x[ii - 1][jj][kk]) * lat.inv_d[0];
                GradNZ = (nd_x[ii][jj][kk] - nd_x[ii][jj][kk - 1]) * lat.inv_d[2];

                VdotGradN = Vx1 * GradNX + Vy1 * Grad_Div + Vz1 * GradNZ;
              }
              VdotGradN = 0.0;

              if (!fdtd::materials::model_spec().harmonic) {
                jy_extra[mmy][ii][jj][kk] = fdtd::materials::electron_gas(mmy).nl4 * (1.0 / (E0 * (TOTND1))) * (VdotGrad + Vy1 * (DivJ -
                    (VdotGradN + VdotGradNeq) / (TOTND1))) - fdtd::materials::electron_gas(mmy).nl4 * (MU0 * E0 / effective_electron_mass[mmy]) * VcrossH +
                  fdtd::materials::electron_gas(mmy).nl5 * pow(TOTND1, 2.0 / 3.0) * Grad_Div;
              } else {
                jy_extra[mmy][ii][jj][kk] = fdtd::materials::electron_gas(mmy).nl4 * inverse_chargedens[mmy] * (VdotGrad + Vy1 * (DivJ)) -
                  fdtd::materials::electron_gas(mmy).nl4 * (MU0 * E0 / effective_electron_mass[mmy]) * VcrossH;
              }
            } else {

              if (with_pressure) Grad_Div = lat.inv_d[1] * (nd_x[ii][jj + 1][kk] - nd_x[ii][jj][kk]);
              else Grad_Div = 0.0;
              ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii][jj + 1][kk]);
              TOTND1 = ND1 + equilibrium_electron_density[mmy];

              Px_this = jx_prev[mmy][ii][jj][kk];
              Py_this = jy_prev[mmy][ii][jj][kk];
              Pz_this = jz_prev[mmy][ii][jj][kk];
              Px_this2 = ds.c[0].j[mmy][ii][jj][kk];
              Py_this2 = ds.c[1].j[mmy][ii][jj][kk];
              Pz_this2 = ds.c[2].j[mmy][ii][jj][kk];
              Px_thisyp1 = jx_prev[mmy][ii][jj + 1][kk];
              Py_thisyp1 = jy_prev[mmy][ii][jj + 1][kk];
              Pz_thisyp1 = jz_prev[mmy][ii][jj + 1][kk];
              Px_this2yp1 = ds.c[0].j[mmy][ii][jj + 1][kk];
              Py_this2yp1 = ds.c[1].j[mmy][ii][jj + 1][kk];
              Pz_this2yp1 = ds.c[2].j[mmy][ii][jj + 1][kk];

              Py_thisym1 = jy_prev[mmy][ii][jj - 1][kk];
              Py_this2ym1 = ds.c[1].j[mmy][ii][jj - 1][kk];

              Px_thisxm1 = jx_prev[mmy][ii - 1][jj][kk];
              Px_this2xm1 = ds.c[0].j[mmy][ii - 1][jj][kk];
              Px_thisxm1yp1 = jx_prev[mmy][ii - 1][jj + 1][kk];
              Px_this2xm1yp1 = ds.c[0].j[mmy][ii - 1][jj + 1][kk];

              Pz_thiszm1 = jz_prev[mmy][ii][jj][kk - 1];
              Pz_this2zm1 = ds.c[2].j[mmy][ii][jj][kk - 1];
              Pz_thiszm1yp1 = jz_prev[mmy][ii][jj + 1][kk - 1];
              Pz_this2zm1yp1 = ds.c[2].j[mmy][ii][jj + 1][kk - 1];

              Vy1 = 0.5 * (Py_this + Py_this2);
              Vx1 = 0.25 * (Px_this + Px_thisyp1 + Px_thisxm1 + Px_thisxm1yp1);
              Vx1 = Vx1 + 0.25 * (Px_this2 + Px_this2yp1 + Px_this2xm1 + Px_this2xm1yp1);
              Vx1 /= 2.0;

              Vz1 = 0.25 * (Pz_this + Pz_thisyp1 + Pz_thiszm1 + Pz_thiszm1yp1);
              Vz1 = Vz1 + 0.25 * (Pz_this2 + Pz_this2yp1 + Pz_this2zm1 + Pz_this2zm1yp1);
              Vz1 /= 2.0;

              Hx1 = 0.5 * (Hx[i][j][k] + Hx[i][j][k - 1]);
              Hz1 = 0.5 * (Hz[i][j][k] + Hz[i - 1][j][k]);
              Hx2 = 0.5 * (hx_prev[i][j][k] + hx_prev[i][j][k - 1]);
              Hz2 = 0.5 * (hz_prev[i][j][k] + hz_prev[i - 1][j][k]);

              Hx1 = 0.5 * (Hx1 + Hx2);
              Hz1 = 0.5 * (Hz1 + Hz2);

              if (with_mag_field) {
                VcrossH = Vz1 * Hx1 - Vx1 * Hz1;
              } else VcrossH = 0.0;

              if (with_convection) {
                
                VdotGrad = 0.25 * (Vx1 * (jy_prev[mmy][ii + 1][jj][kk] + ds.c[1].j[mmy][ii + 1][jj][kk] - jy_prev[mmy][ii - 1][jj][kk] - ds.c[1].j[mmy][ii - 1][jj][kk]) * lat.inv_d[0] +
                  Vy1 * (Py_thisyp1 + Py_this2yp1 - Py_thisym1 - Py_this2ym1) * lat.inv_d[1] +
                  Vz1 * (jy_prev[mmy][ii][jj][kk + 1] + ds.c[1].j[mmy][ii][jj][kk + 1] - jy_prev[mmy][ii][jj][kk - 1] - ds.c[1].j[mmy][ii][jj][kk - 1]) * lat.inv_d[2]);
              } else VdotGrad = 0.0;
              
              DivJ = 0.25 * (lat.inv_d[1] * (Py_thisyp1 + Py_this2yp1 - Py_thisym1 - Py_this2ym1) +
                lat.inv_d[0] * (Px_this + Px_thisyp1 + Px_this2 + Px_this2yp1 -
                  Px_thisxm1 - Px_thisxm1yp1 - Px_this2xm1 - Px_this2xm1yp1) +
                lat.inv_d[2] * (Pz_this + Pz_thisyp1 + Pz_this2 + Pz_this2yp1 -
                  Pz_thiszm1 - Pz_thiszm1yp1 - Pz_this2zm1 - Pz_this2zm1yp1));

              VdotGradNeq = 0.0;

              if (!fdtd::materials::model_spec().harmonic) {
                GradNX = (nd_x[ii + 1][jj + 1][kk] + nd_x[ii + 1][jj][kk] - nd_x[ii - 1][jj + 1][kk] - nd_x[ii - 1][jj][kk]) * 0.25 * lat.inv_d[0];
                GradNZ = (nd_x[ii][jj + 1][kk + 1] + nd_x[ii][jj][kk + 1] - nd_x[ii][jj + 1][kk - 1] - nd_x[ii][jj][kk - 1]) * 0.25 * lat.inv_d[2];

                VdotGradN = Vx1 * GradNX + Vy1 * Grad_Div + Vz1 * GradNZ;
              }

              if (!fdtd::materials::model_spec().harmonic) {
                jy_extra[mmy][ii][jj][kk] = fdtd::materials::electron_gas(mmy).nl4 * (1.0 / (E0 * (TOTND1))) * (VdotGrad + Vy1 * (DivJ -
                  (VdotGradN + VdotGradNeq) / (TOTND1))) - fdtd::materials::electron_gas(mmy).nl4 * (MU0 * E0 / effective_electron_mass[mmy]) * VcrossH; //+ fdtd::materials::electron_gas(mmy).nl5*pow(TOTND1,2.0/3.0)*Grad_Div;

                // Py_NL3[ii][jj][kk] = fdtd::materials::electron_gas(mmy).nl4*inverse_chargedens[mmy] * (VdotGrad + Vy1*DivJ)
              } else {
                jy_extra[mmy][ii][jj][kk] = fdtd::materials::electron_gas(mmy).nl4 * inverse_chargedens[mmy] * (VdotGrad + Vy1 * (DivJ)) -
                  fdtd::materials::electron_gas(mmy).nl4 * (MU0 * E0 / effective_electron_mass[mmy]) * VcrossH;
              }
            }
          }
          
          else if (mmy >= NonLinFreeElectron && mmy < NonlinearDrudeLorentz) {
            if (with_pressure) Grad_Div = lat.inv_d[1] * (nd_x[ii][jj + 1][kk] - nd_x[ii][jj][kk]);
            else Grad_Div = 0.0;
            ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii][jj + 1][kk]);
            TOTND1 = ND1 + equilibrium_electron_density[mmy];

            Px_this = jx_prev[mmy][ii][jj][kk];
            Py_this = jy_prev[mmy][ii][jj][kk];
            Pz_this = jz_prev[mmy][ii][jj][kk];
            Px_this2 = ds.c[0].j[mmy][ii][jj][kk];
            Py_this2 = ds.c[1].j[mmy][ii][jj][kk];
            Pz_this2 = ds.c[2].j[mmy][ii][jj][kk];
            Px_thisyp1 = jx_prev[mmy][ii][jj + 1][kk];
            Py_thisyp1 = jy_prev[mmy][ii][jj + 1][kk];
            Pz_thisyp1 = jz_prev[mmy][ii][jj + 1][kk];
            Px_this2yp1 = ds.c[0].j[mmy][ii][jj + 1][kk];
            Py_this2yp1 = ds.c[1].j[mmy][ii][jj + 1][kk];
            Pz_this2yp1 = ds.c[2].j[mmy][ii][jj + 1][kk];

            Py_thisym1 = jy_prev[mmy][ii][jj - 1][kk];
            Py_this2ym1 = ds.c[1].j[mmy][ii][jj - 1][kk];

            Px_thisxm1 = jx_prev[mmy][ii - 1][jj][kk];
            Px_this2xm1 = ds.c[0].j[mmy][ii - 1][jj][kk];
            Px_thisxm1yp1 = jx_prev[mmy][ii - 1][jj + 1][kk];
            Px_this2xm1yp1 = ds.c[0].j[mmy][ii - 1][jj + 1][kk];

            Pz_thiszm1 = jz_prev[mmy][ii][jj][kk - 1];
            Pz_this2zm1 = ds.c[2].j[mmy][ii][jj][kk - 1];
            Pz_thiszm1yp1 = jz_prev[mmy][ii][jj + 1][kk - 1];
            Pz_this2zm1yp1 = ds.c[2].j[mmy][ii][jj + 1][kk - 1];

            Vy1 = 0.5 * (Py_this + Py_this2);
            Vx1 = 0.25 * (Px_this + Px_thisyp1 + Px_thisxm1 + Px_thisxm1yp1);
            Vx1 = Vx1 + 0.25 * (Px_this2 + Px_this2yp1 + Px_this2xm1 + Px_this2xm1yp1);
            Vx1 /= 2.0;

            Vz1 = 0.25 * (Pz_this + Pz_thisyp1 + Pz_thiszm1 + Pz_thiszm1yp1);
            Vz1 = Vz1 + 0.25 * (Pz_this2 + Pz_this2yp1 + Pz_this2zm1 + Pz_this2zm1yp1);
            Vz1 /= 2.0;

            Hx1 = 0.5 * (Hx[i][j][k] + Hx[i][j][k - 1]);
            Hz1 = 0.5 * (Hz[i][j][k] + Hz[i - 1][j][k]);
            Hx2 = 0.5 * (hx_prev[i][j][k] + hx_prev[i][j][k - 1]);
            Hz2 = 0.5 * (hz_prev[i][j][k] + hz_prev[i - 1][j][k]);

            Hx1 = 0.5 * (Hx1 + Hx2);
            Hz1 = 0.5 * (Hz1 + Hz2);

            if (with_mag_field) {
              VcrossH = Vz1 * Hx1 - Vx1 * Hz1;
            } else VcrossH = 0.0;

            if (with_convection) {

              VdotGrad = 0.25 * (Vx1 * (jy_prev[mmy][ii + 1][jj][kk] + ds.c[1].j[mmy][ii + 1][jj][kk] - jy_prev[mmy][ii - 1][jj][kk] - ds.c[1].j[mmy][ii - 1][jj][kk]) * lat.inv_d[0] +
                Vy1 * (Py_thisyp1 + Py_this2yp1 - Py_thisym1 - Py_this2ym1) * lat.inv_d[1] +
                Vz1 * (jy_prev[mmy][ii][jj][kk + 1] + ds.c[1].j[mmy][ii][jj][kk + 1] - jy_prev[mmy][ii][jj][kk - 1] - ds.c[1].j[mmy][ii][jj][kk - 1]) * lat.inv_d[2]);
            } else VdotGrad = 0.0;

            DivJ = 0.25 * (lat.inv_d[1] * (Py_thisyp1 + Py_this2yp1 - Py_thisym1 - Py_this2ym1) +
              lat.inv_d[0] * (Px_this + Px_thisyp1 + Px_this2 + Px_this2yp1 -
                Px_thisxm1 - Px_thisxm1yp1 - Px_this2xm1 - Px_this2xm1yp1) +
              lat.inv_d[2] * (Pz_this + Pz_thisyp1 + Pz_this2 + Pz_this2yp1 -
                Pz_thiszm1 - Pz_thiszm1yp1 - Pz_this2zm1 - Pz_this2zm1yp1));

            VdotGradNeq = 0.0;
            VdotGradN = 0.0;

            NONLINDrude4TTM = 1.0 / (1.0 / fdtd::materials::model_spec().dt + 0.5 * fdtd::materials::thermal::damping_constant(fdtd::materials::thermal::electron_temp[i][j][k], fdtd::materials::thermal::lattice_temp[i][j][k]));

            jy_extra[mmy][ii][jj][kk] = NONLINDrude4TTM * (1.0 / (E0 * (TOTND1))) * (VdotGrad + Vy1 * (DivJ -
              (VdotGradN + VdotGradNeq) / (TOTND1))) - NONLINDrude4TTM * ((MU0 * EP0 * fdtd::materials::thermal::plasma_frequency(fdtd::materials::thermal::electron_temp[i][j][k]) * fdtd::materials::thermal::plasma_frequency(fdtd::materials::thermal::electron_temp[i][j][k])) / (E0 * TOTND1)) * VcrossH; //+ fdtd::materials::electron_gas(mmy).nl5*pow(TOTND1,2.0/3.0)*Grad_Div;

          }
        }

      }

    }
  }

  return;
}

/**
 * @brief 1st step of update algorithm for hydrodynamic model
 *
 * @details Updates jz using a "linear/local" approximation - the first step in
 * the hydrodyanmic update procedure. Also supports the nonlinear Drude-Lorentz
 * model + TTM.
 */
void fdtd::materials::hydro::update_pz()
{
	auto &ds = fdtd::materials::dispersive_state();
	int i, j, k, u, v, w, U, V, W, mz, mmz; //
	int ii, jj, kk;
	double Grad_Div, Grad_Div2;
	double ND1, powerapprox;
	double gamma;



	for (i = 1; i <= subdomain.owned_cells_x; ++i)
	{
		ii = i + 1;

		if (i >= subdomain.i_cpml_min && i <= subdomain.i_cpml_max)
		{

			U = 1;
			u = i - subdomain.inner_org_x;
		}
		else
		{

			if (i < subdomain.i_cpml_min)
			{

                U = 0;
			}
			else
			{

                U = 2;
			}

		}

		for (j = 1; j <= subdomain.owned_cells_y; ++j)
		{
			jj = j + 1;

			if (j >= subdomain.j_cpml_min && j <= subdomain.j_cpml_max)
			{

				V = 1;
				v = j - subdomain.inner_org_y;
			}
			else
			{

				if (j < subdomain.j_cpml_min)
				{

                    V = 0;
				}
				else
				{

                    V = 2;
				}

			}

            for (k = 0; k < subdomain.cells_z; ++k) {
                kk = k + 1;

				if (k >= subdomain.k_cpml_min && k < subdomain.k_cpml_max)
				{

					W = 1;
					w = k - subdomain.inner_org_z;
				}
				else
				{

					if (k < subdomain.k_cpml_min)
					{

						W = 0;
					}
					else
					{

						W = 2;
					}
				}

				if ((U == 1) && (V == 1) && (W == 1))
				{
					mz = int(matind_z[u][v][w]);
				}
				else
				{

					if (V == 2)
						mz = fdtd::materials::model_spec().medium_second;

					else if (V == 0)
						mz = fdtd::materials::model_spec().medium_first;

					else
						mz = background_y[j + subdomain.j_for_uni];
				}

				// if(mz==nl_code){mz=NL_mat; mz_=nl_code;} //(SHG || THG) &&
				if (mz >= DispMateOffsetD && mz < NRMAXMATERIALS)
				{
					mmz = mz - DispMateOffsetD;

					if (mmz >= ModLorentz && mmz < NonLinFreeElectron)
					{ // Free Electron Hydrodynamics model
						if (with_pressure)
							Grad_Div = lat.inv_d[2] * (nd_x[ii][jj][kk + 1] - nd_x[ii][jj][kk]);
						else
							Grad_Div = 0.0;
						ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii][jj][kk + 1]);
						if (!fdtd::materials::model_spec().harmonic)
						{
							ds.c[2].j[mmz][ii][jj][kk] = fdtd::materials::electron_gas(mmz).nl1 * jz_prev[mmz][ii][jj][kk] + (fdtd::materials::electron_gas(mmz).nl2 + fdtd::materials::electron_gas(mmz).nl3 * ND1) * Ez[i][j][k];
						}
						else
						{
							if (with_pressure)
								Grad_Div2 = lat.inv_d[2] * (nd_y[ii][jj][kk + 1] - nd_y[ii][jj][kk]);
							else
								Grad_Div2 = 0.0;
							powerapprox = powerterm[mmz] * (1.0 + (2.0 / 3.0) * ND1 * inverse_edens[mmz]);

							ds.c[2].j[mmz][ii][jj][kk] = fdtd::materials::electron_gas(mmz).nl1 * jz_prev[mmz][ii][jj][kk] + (fdtd::materials::electron_gas(mmz).nl2) * Ez[i][j][k] + fdtd::materials::electron_gas(mmz).nl5 * powerterm[mmz] * Grad_Div;
							ds.c[2].j_nl[mmz][ii][jj][kk] = fdtd::materials::electron_gas(mmz).nl1 * ds.c[2].j_nl[mmz][ii][jj][kk] + (fdtd::materials::electron_gas(mmz).nl2 * ez_nl[i][j][k] + fdtd::materials::electron_gas(mmz).nl3 * ND1 * Ez[i][j][k]) + fdtd::materials::electron_gas(mmz).nl5 * powerterm[mmz] * Grad_Div2 + fdtd::materials::electron_gas(mmz).nl5 * powerapprox * Grad_Div;
						}
					}

					else if (mmz >= NonLinFreeElectron && mmz < NonlinearDrudeLorentz)
					{

						gamma = fdtd::materials::thermal::damping_constant(fdtd::materials::thermal::electron_temp[i][j][k], fdtd::materials::thermal::lattice_temp[i][j][k]);
						ds.c[2].j[mmz][ii][jj][kk] = ((1.0 / fdtd::materials::model_spec().dt - 0.5 * gamma) / (1.0 / fdtd::materials::model_spec().dt + 0.5 * gamma)) * jz_prev[mmz][ii][jj][kk] + (EP0 * fdtd::materials::thermal::plasma_frequency((fdtd::materials::thermal::electron_temp[i][j][k])) * fdtd::materials::thermal::plasma_frequency((fdtd::materials::thermal::electron_temp[i][j][k])) / (1.0 / fdtd::materials::model_spec().dt + 0.5 * gamma)) * Ez[i][j][k];
					}
				}
			}
		}
	}

	return;
}

/**
 * @brief 2nd step of update algorithm for hydrodynamic model
 *
 * @details Calculates the nonlinear correction to jz, see docs
 */
void fdtd::materials::hydro::update_pz_2() {
  auto &ds = fdtd::materials::dispersive_state();
  int i, j, k, u, v, w, U, V, W, mz, mmz; //
  int ii, jj, kk;
  double Grad_Div;
  double Vx1, Vy1, Vz1, Hx1, Hy1, Hx2, Hy2, VdotGrad, VcrossH, ND1, DivJ, VdotGradN, VdotGradNeq;
  double GradNX, GradNY, TOTND1;
  double Px_this, Py_this, Pz_this, Px_this2, Py_this2, Pz_this2, Px_thiszp1, Py_thiszp1, Pz_thiszp1, Px_this2zp1, Py_this2zp1, Pz_this2zp1, Pz_thiszm1, Pz_this2zm1, Px_thisxm1, Px_this2xm1, Px_thisxm1zp1, Px_this2xm1zp1, Py_thisym1, Py_this2ym1, Py_thisym1zp1, Py_this2ym1zp1;

  double NONLINDrude4TTM;

  for (i = 1; i <= subdomain.owned_cells_x; ++i) {
    ii = i + 1;

    if (i >= subdomain.i_cpml_min && i <= subdomain.i_cpml_max) {

      U = 1;
      u = i - subdomain.inner_org_x;

    } else {

      if (i < subdomain.i_cpml_min) {

        U = 0;

      } else {

        U = 2;
      }

    }

    for (j = 1; j <= subdomain.owned_cells_y; ++j) {
      jj = j + 1;

      if (j >= subdomain.j_cpml_min && j <= subdomain.j_cpml_max) {

        V = 1;
        v = j - subdomain.inner_org_y;

      } else {

        if (j < subdomain.j_cpml_min) {

          V = 0;

        } else {

          V = 2;
        }

      }

      for (k = 0; k < subdomain.cells_z; ++k) {
        kk = k + 1;

        if (k >= subdomain.k_cpml_min && k < subdomain.k_cpml_max) {

          W = 1;
          w = k - subdomain.inner_org_z;

        } else {

          if (k < subdomain.k_cpml_min) {

            W = 0;

          } else {

            W = 2;
          }
        }

        if ((U == 1) && (V == 1) && (W == 1)) {
          mz = int(matind_z[u][v][w]);
        } else {

          if (V == 2) mz = fdtd::materials::model_spec().medium_second;

          else if (V == 0) mz = fdtd::materials::model_spec().medium_first;

          else mz = background_y[j + subdomain.j_for_uni];

        }
        //if(mz==nl_code){mz=NL_mat; mz_=nl_code;} //(SHG || THG) &&
        if (mz >= DispMateOffsetD && mz < NRMAXMATERIALS) {
          mmz = mz - DispMateOffsetD;

          if (mmz >= ModLorentz && mmz < NonLinFreeElectron) { //Free Electron Hydrodynamics model

            if (use_backwards_differences) {
              if (with_pressure) Grad_Div = lat.inv_d[2] * (nd_x[ii][jj][kk + 1] - nd_x[ii][jj][kk]);
              else Grad_Div = 0.0;
              ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii][jj][kk + 1]);
              TOTND1 = ND1 + equilibrium_electron_density[mmz];

              Px_this = jx_prev[mmz][ii][jj][kk];
              Py_this = jy_prev[mmz][ii][jj][kk];
              Pz_this = jz_prev[mmz][ii][jj][kk];

              Px_thiszp1 = jx_prev[mmz][ii][jj][kk + 1];
              Py_thiszp1 = jy_prev[mmz][ii][jj][kk + 1];
              Pz_thiszp1 = jz_prev[mmz][ii][jj][kk + 1];

              Pz_thiszm1 = jz_prev[mmz][ii][jj][kk - 1];

              Px_thisxm1 = jx_prev[mmz][ii - 1][jj][kk];
              Px_thisxm1zp1 = jx_prev[mmz][ii - 1][jj][kk + 1];

              Py_thisym1 = jy_prev[mmz][ii][jj - 1][kk];
              Py_thisym1zp1 = jy_prev[mmz][ii][jj - 1][kk + 1];

              Vz1 = (Pz_this);
              Vx1 = (Px_this);
              Vy1 = (Py_this);

              Hy1 = (Hy[i][j][k]);
              Hx1 = (Hx[i][j][k]);
              Hy2 = (hy_prev[i][j][k]);
              Hx2 = (hx_prev[i][j][k]);

              Hx1 = 0.5 * (Hx1 + Hx2);
              Hy1 = 0.5 * (Hy1 + Hy2);

              if (with_mag_field) {
                VcrossH = Vx1 * Hy1 - Vy1 * Hx1;
              } else VcrossH = 0.0;

              if (with_convection) {
                
                VdotGrad = (Vx1 * (Pz_this - jz_prev[mmz][ii - 1][jj][kk]) * lat.inv_d[0] +
                  Vy1 * (Pz_this - jz_prev[mmz][ii][jj - 1][kk]) * lat.inv_d[1] +
                  Vz1 * (Pz_this - Pz_thiszm1) * lat.inv_d[2]);
              } else VdotGrad = 0.0;

              DivJ = (lat.inv_d[2] * (Pz_this - Pz_thiszm1) +
                lat.inv_d[0] * (Px_this - Px_thisxm1) +
                lat.inv_d[1] * (Py_this - Py_thisym1));
              
              VdotGradNeq = 0.0;

              if (!fdtd::materials::model_spec().harmonic) {
                GradNX = (nd_x[ii][jj][kk] - nd_x[ii - 1][jj][kk]) * lat.inv_d[0];
                GradNY = (nd_x[ii][jj][kk] - nd_x[ii][jj - 1][kk]) * lat.inv_d[1];
                VdotGradN = Vx1 * GradNX + Vy1 * GradNY + Vz1 * Grad_Div;
              }
              VdotGradN = 0.0;

              if (!fdtd::materials::model_spec().harmonic) {
                jz_extra[mmz][ii][jj][kk] = fdtd::materials::electron_gas(mmz).nl4 * (1.0 / (E0 * (TOTND1))) * (VdotGrad + Vz1 * (DivJ -
                    (VdotGradN + VdotGradNeq) / (TOTND1))) - fdtd::materials::electron_gas(mmz).nl4 * (MU0 * E0 / effective_electron_mass[mmz]) * VcrossH +
                  fdtd::materials::electron_gas(mmz).nl5 * pow(TOTND1, 2.0 / 3.0) * Grad_Div;
              } else {
                jz_extra[mmz][ii][jj][kk] = fdtd::materials::electron_gas(mmz).nl4 * inverse_chargedens[mmz] * (VdotGrad + Vz1 * (DivJ)) -
                  fdtd::materials::electron_gas(mmz).nl4 * (MU0 * E0 / effective_electron_mass[mmz]) * VcrossH;
              }
            } else {
              if (with_pressure) Grad_Div = lat.inv_d[2] * (nd_x[ii][jj][kk + 1] - nd_x[ii][jj][kk]);
              else Grad_Div = 0.0;
              ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii][jj][kk + 1]);
              TOTND1 = ND1 + equilibrium_electron_density[mmz];

              Px_this = jz_prev[mmz][ii][jj][kk];
              Py_this = jy_prev[mmz][ii][jj][kk];
              Pz_this = jz_prev[mmz][ii][jj][kk];
              Px_this2 = ds.c[0].j[mmz][ii][jj][kk];
              Py_this2 = ds.c[1].j[mmz][ii][jj][kk];
              Pz_this2 = ds.c[2].j[mmz][ii][jj][kk];
              Px_thiszp1 = jx_prev[mmz][ii][jj][kk + 1];
              Py_thiszp1 = jy_prev[mmz][ii][jj][kk + 1];
              Pz_thiszp1 = jz_prev[mmz][ii][jj][kk + 1];
              Px_this2zp1 = ds.c[0].j[mmz][ii][jj][kk + 1];
              Py_this2zp1 = ds.c[1].j[mmz][ii][jj][kk + 1];
              Pz_this2zp1 = ds.c[2].j[mmz][ii][jj][kk + 1];

              Pz_thiszm1 = jz_prev[mmz][ii][jj][kk - 1];
              Pz_this2zm1 = ds.c[2].j[mmz][ii][jj][kk - 1];

              Px_thisxm1 = jx_prev[mmz][ii - 1][jj][kk];
              Px_this2xm1 = ds.c[0].j[mmz][ii - 1][jj][kk];
              Px_thisxm1zp1 = jx_prev[mmz][ii - 1][jj][kk + 1];
              Px_this2xm1zp1 = ds.c[0].j[mmz][ii - 1][jj][kk + 1];

              Py_thisym1 = jy_prev[mmz][ii][jj - 1][kk];
              Py_this2ym1 = ds.c[1].j[mmz][ii][jj - 1][kk];
              Py_thisym1zp1 = jy_prev[mmz][ii][jj - 1][kk + 1];
              Py_this2ym1zp1 = ds.c[1].j[mmz][ii][jj - 1][kk + 1];

              Vz1 = 0.5 * (Pz_this + Pz_this2);
              Vx1 = 0.25 * (Px_this + Px_thiszp1 + Px_thisxm1zp1 + Px_thisxm1);
              Vx1 = Vx1 + 0.25 * (Px_this2 + Px_this2zp1 + Px_this2xm1zp1 + Px_this2xm1);
              Vx1 /= 2.0;
              Vy1 = 0.25 * (Py_this + Py_thiszp1 + Py_thisym1zp1 + Py_thisym1);
              Vy1 = Vy1 + 0.25 * (Py_this2 + Py_this2zp1 + Py_this2ym1zp1 + Py_this2ym1);
              Vy1 /= 2.0;

              Hy1 = 0.5 * (Hy[i - 1][j][k] + Hy[i][j][k]);
              Hx1 = 0.5 * (Hx[i][j][k] + Hx[i][j - 1][k]);
              Hy2 = 0.5 * (hy_prev[i - 1][j][k] + hy_prev[i][j][k]);
              Hx2 = 0.5 * (hx_prev[i][j][k] + hx_prev[i][j - 1][k]);

              Hx1 = 0.5 * (Hx1 + Hx2);
              Hy1 = 0.5 * (Hy1 + Hy2);

              if (with_mag_field) {
                VcrossH = Vx1 * Hy1 - Vy1 * Hx1;
              } else VcrossH = 0.0;

              if (with_convection) {
               
                VdotGrad = 0.25 * (Vx1 * (jz_prev[mmz][ii + 1][jj][kk] + ds.c[2].j[mmz][ii + 1][jj][kk] - jz_prev[mmz][ii - 1][jj][kk] - ds.c[2].j[mmz][ii - 1][jj][kk]) * lat.inv_d[0] +
                  Vy1 * (jz_prev[mmz][ii][jj + 1][kk] + ds.c[2].j[mmz][ii][jj + 1][kk] - jz_prev[mmz][ii][jj - 1][kk] - ds.c[2].j[mmz][ii][jj - 1][kk]) * lat.inv_d[1] +
                  Vz1 * (Pz_thiszp1 + Pz_this2zp1 - Pz_thiszm1 - Pz_this2zm1) * lat.inv_d[2]);
              } else VdotGrad = 0.0;

              DivJ = 0.25 * (lat.inv_d[2] * (Pz_thiszp1 + Pz_this2zp1 - Pz_thiszm1 - Pz_this2zp1) +
                lat.inv_d[0] * (Px_this + Px_thiszp1 + Px_this2 + Px_this2zp1 -
                  Px_thisxm1 - Px_thisxm1zp1 - Px_this2xm1 - Px_this2xm1zp1) +
                lat.inv_d[1] * (Py_this + Py_thiszp1 + Py_this2 + Py_this2zp1 -
                  Py_thisym1 - Py_thisym1zp1 - Py_this2ym1 - Py_this2ym1zp1));
              
              VdotGradNeq = 0.0;

              if (!fdtd::materials::model_spec().harmonic) {
                GradNX = (nd_x[ii + 1][jj][kk + 1] + nd_x[ii + 1][jj][kk] - nd_x[ii - 1][jj][kk + 1] - nd_x[ii - 1][jj][kk]) * 0.25 * lat.inv_d[0];
                GradNY = (nd_x[ii][jj + 1][kk + 1] + nd_x[ii][jj + 1][kk] - nd_x[ii][jj - 1][kk + 1] - nd_x[ii][jj - 1][kk]) * 0.25 * lat.inv_d[1];
                VdotGradN = Vx1 * GradNX + Vy1 * GradNY + Vz1 * Grad_Div;
              }

              if (!fdtd::materials::model_spec().harmonic) {
                jz_extra[mmz][ii][jj][kk] = fdtd::materials::electron_gas(mmz).nl4 * (1.0 / (E0 * (TOTND1))) * (VdotGrad + Vz1 * (DivJ -
                  (VdotGradN + VdotGradNeq) / (TOTND1))) - fdtd::materials::electron_gas(mmz).nl4 * (MU0 * E0 / effective_electron_mass[mmz]) * VcrossH; 

              } else {
                jz_extra[mmz][ii][jj][kk] = fdtd::materials::electron_gas(mmz).nl4 * inverse_chargedens[mmz] * (VdotGrad + Vz1 * (DivJ)) -
                  fdtd::materials::electron_gas(mmz).nl4 * (MU0 * E0 / effective_electron_mass[mmz]) * VcrossH;
              }
            }

          }

          else if (mmz >= NonLinFreeElectron && mmz < NonlinearDrudeLorentz) {
            if (with_pressure) Grad_Div = lat.inv_d[2] * (nd_x[ii][jj][kk + 1] - nd_x[ii][jj][kk]);
            else Grad_Div = 0.0;
            ND1 = 0.5 * (nd_x[ii][jj][kk] + nd_x[ii][jj][kk + 1]);
            TOTND1 = ND1 + equilibrium_electron_density[mmz];

            Px_this = jz_prev[mmz][ii][jj][kk];
            Py_this = jy_prev[mmz][ii][jj][kk];
            Pz_this = jz_prev[mmz][ii][jj][kk];
            Px_this2 = ds.c[0].j[mmz][ii][jj][kk];
            Py_this2 = ds.c[1].j[mmz][ii][jj][kk];
            Pz_this2 = ds.c[2].j[mmz][ii][jj][kk];
            Px_thiszp1 = jx_prev[mmz][ii][jj][kk + 1];
            Py_thiszp1 = jy_prev[mmz][ii][jj][kk + 1];
            Pz_thiszp1 = jz_prev[mmz][ii][jj][kk + 1];
            Px_this2zp1 = ds.c[0].j[mmz][ii][jj][kk + 1];
            Py_this2zp1 = ds.c[1].j[mmz][ii][jj][kk + 1];
            Pz_this2zp1 = ds.c[2].j[mmz][ii][jj][kk + 1];

            Pz_thiszm1 = jz_prev[mmz][ii][jj][kk - 1];
            Pz_this2zm1 = ds.c[2].j[mmz][ii][jj][kk - 1];

            Px_thisxm1 = jx_prev[mmz][ii - 1][jj][kk];
            Px_this2xm1 = ds.c[0].j[mmz][ii - 1][jj][kk];
            Px_thisxm1zp1 = jx_prev[mmz][ii - 1][jj][kk + 1];
            Px_this2xm1zp1 = ds.c[0].j[mmz][ii - 1][jj][kk + 1];

            Py_thisym1 = jy_prev[mmz][ii][jj - 1][kk];
            Py_this2ym1 = ds.c[1].j[mmz][ii][jj - 1][kk];
            Py_thisym1zp1 = jy_prev[mmz][ii][jj - 1][kk + 1];
            Py_this2ym1zp1 = ds.c[1].j[mmz][ii][jj - 1][kk + 1];

            Vz1 = 0.5 * (Pz_this + Pz_this2);
            Vx1 = 0.25 * (Px_this + Px_thiszp1 + Px_thisxm1zp1 + Px_thisxm1);
            Vx1 = Vx1 + 0.25 * (Px_this2 + Px_this2zp1 + Px_this2xm1zp1 + Px_this2xm1);
            Vx1 /= 2.0;
            Vy1 = 0.25 * (Py_this + Py_thiszp1 + Py_thisym1zp1 + Py_thisym1);
            Vy1 = Vy1 + 0.25 * (Py_this2 + Py_this2zp1 + Py_this2ym1zp1 + Py_this2ym1);
            Vy1 /= 2.0;

            Hy1 = 0.5 * (Hy[i - 1][j][k] + Hy[i][j][k]);
            Hx1 = 0.5 * (Hx[i][j][k] + Hx[i][j - 1][k]);
            Hy2 = 0.5 * (hy_prev[i - 1][j][k] + hy_prev[i][j][k]);
            Hx2 = 0.5 * (hx_prev[i][j][k] + hx_prev[i][j - 1][k]);

            Hx1 = 0.5 * (Hx1 + Hx2);
            Hy1 = 0.5 * (Hy1 + Hy2);

            if (with_mag_field) {
              VcrossH = Vx1 * Hy1 - Vy1 * Hx1;
            } else VcrossH = 0.0;

            if (with_convection) {
              VdotGrad = 0.25 * (Vx1 * (jz_prev[mmz][ii + 1][jj][kk] + ds.c[2].j[mmz][ii + 1][jj][kk] - jz_prev[mmz][ii - 1][jj][kk] - ds.c[2].j[mmz][ii - 1][jj][kk]) * lat.inv_d[0] +
                Vy1 * (jz_prev[mmz][ii][jj + 1][kk] + ds.c[2].j[mmz][ii][jj + 1][kk] - jz_prev[mmz][ii][jj - 1][kk] - ds.c[2].j[mmz][ii][jj - 1][kk]) * lat.inv_d[1] +
                Vz1 * (Pz_thiszp1 + Pz_this2zp1 - Pz_thiszm1 - Pz_this2zm1) * lat.inv_d[2]);
            } else VdotGrad = 0.0;

            DivJ = 0.25 * (lat.inv_d[2] * (Pz_thiszp1 + Pz_this2zp1 - Pz_thiszm1 - Pz_this2zp1) +
              lat.inv_d[0] * (Px_this + Px_thiszp1 + Px_this2 + Px_this2zp1 -
                Px_thisxm1 - Px_thisxm1zp1 - Px_this2xm1 - Px_this2xm1zp1) +
              lat.inv_d[1] * (Py_this + Py_thiszp1 + Py_this2 + Py_this2zp1 -
                Py_thisym1 - Py_thisym1zp1 - Py_this2ym1 - Py_this2ym1zp1));

            VdotGradNeq = 0.0;
            VdotGradN = 0.0;
            NONLINDrude4TTM = 1.0 / (1.0 / fdtd::materials::model_spec().dt + 0.5 * fdtd::materials::thermal::damping_constant(fdtd::materials::thermal::electron_temp[i][j][k], fdtd::materials::thermal::lattice_temp[i][j][k]));

            jz_extra[mmz][ii][jj][kk] = NONLINDrude4TTM * (1.0 / (E0 * (TOTND1))) * (VdotGrad + Vz1 * (DivJ -
              (VdotGradN + VdotGradNeq) / (TOTND1))) - NONLINDrude4TTM * ((MU0 * EP0 * fdtd::materials::thermal::plasma_frequency(fdtd::materials::thermal::electron_temp[i][j][k]) * fdtd::materials::thermal::plasma_frequency(fdtd::materials::thermal::electron_temp[i][j][k])) / (E0 * TOTND1)) * VcrossH;

          }

        }

      }
    }
  }

  return;
}
