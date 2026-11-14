#include "fields/fields.h"
#include "io/text.h"
#include "materials/hydro/state.h"
#include "materials/hydro/hydro.h"
#include "grid/halo_exchange.h"
#include "grid/subdomain.h"
#include "grid/lattice.h"
#include "materials/spec.h"
using fdtd::constants::E0;
using fdtd::constants::EP0;
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::grid::subdomain;

namespace
{
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::materials::hydro; // transitional: module state, hydro/state.h

void fdtd::materials::hydro::update_nonlocal() {
	int i, j, k;


	fdtd::grid::exchange_margins(px_nl, py_nl, pz_nl);

	if (!fdtd::materials::model_spec().hydrodynamics) {
		for (i = 0; i <= subdomain.cells_x + 2; i++) {
			for (j = 0; j <= subdomain.cells_y + 2; j++) {
				for (k = 0; k <= subdomain.cells_z + 2; k++) {

					px_nl_2[i][j][k] = px_nl[i][j][k];
					py_nl_2[i][j][k] = py_nl[i][j][k];
					pz_nl_2[i][j][k] = pz_nl[i][j][k];
					if (fdtd::materials::model_spec().gnor && i <= subdomain.cells_x && j <= subdomain.cells_y && k <= subdomain.cells_z) {
						grad_div_x_n_1[i][j][k] = grad_div_x_n[i][j][k];
						grad_div_x_n[i][j][k] = grad_div_x[i][j][k];
						grad_div_y_n_1[i][j][k] = grad_div_y_n[i][j][k];
						grad_div_y_n[i][j][k] = grad_div_y[i][j][k];
						grad_div_z_n_1[i][j][k] = grad_div_z_n[i][j][k];
						grad_div_z_n[i][j][k] = grad_div_z[i][j][k];
					}
				}
			}
		}
	}
	return;
}

void fdtd::materials::hydro::update_electron_density_local() {
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
			}
		}
	}

	fdtd::grid::exchange_margins(ex_2, ey_2, ez_2);
	fdtd::grid::exchange_margins(px_cp_2, py_cp_2, pz_cp_2);

	if (!subdomain.cpml_only) {
		for (i = 0; i <= subdomain.cells_x + 1; i++) {
			ii = i + 1;
			for (j = 0; j <= subdomain.cells_y + 1; j++) {
				jj = j + 1;
				for (k = 0; k <= subdomain.cells_z + 1; k++) {
					kk = k + 1;
					nd_x[ii][jj][kk] = (-1.0 / E0) * (EP0 * (lat.inv_d[0] * (ex_2[ii][jj][kk] * epsinf_x[ii][jj][kk] - ex_2[ii - 1][jj][kk] * epsinf_x[ii - 1][jj][kk]) +
						lat.inv_d[1] * (ey_2[ii][jj][kk] * epsinf_y[ii][jj][kk] - ey_2[ii][jj - 1][kk] * epsinf_y[ii][jj - 1][kk]) +
						lat.inv_d[2] * (ez_2[ii][jj][kk] * epsinf_z[ii][jj][kk] - ez_2[ii][jj][kk - 1] * epsinf_z[ii][jj][kk - 1])) +
						lat.inv_d[0] * (px_cp_2[ii][jj][kk] - px_cp_2[ii - 1][jj][kk]) +
						lat.inv_d[1] * (py_cp_2[ii][jj][kk] - py_cp_2[ii][jj - 1][kk]) +
						lat.inv_d[2] * (pz_cp_2[ii][jj][kk] - pz_cp_2[ii][jj][kk - 1]));
					nd_z[ii][jj][kk] = (-1.0 / E0) * (EP0 * (lat.inv_d[0] * (ex_2[ii][jj][kk] - ex_2[ii - 1][jj][kk]) + lat.inv_d[1] * (ey_2[ii][jj][kk] - ey_2[ii][jj - 1][kk]) + lat.inv_d[2] * (ez_2[ii][jj][kk] - ez_2[ii][jj][kk - 1])));
				}
			}
		}
	}
	fdtd::grid::exchange_margin(nd_x);
	fdtd::grid::exchange_margin(nd_z);
}