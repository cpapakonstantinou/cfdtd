/**
 * @file geometry/composite.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-01-19
 * 
 * @brief Basic shape creation functions
 * 
 * @note must be gone over and organized
 * 
 * @copyright Copyright 2022
 */
#include "materials/material_map.h"
#include "grid/subdomain.h"
#include "sources/signal.h"
#include "sources/aux_grid.h"
#include "grid/state.h"
#include "sources/state.h"
#include "io/text.h"
#include "geometry/state.h"
#include "geometry/builders.h"
#include "grid/lattice.h"
#include "geometry/scene.h"
#include "geometry/shapes.h"
using fdtd::constants::PI;
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;
using fdtd::materials::matind_x;
using fdtd::materials::matind_y;
using fdtd::materials::matind_z;
using fdtd::materials::background_y;

namespace
{
const fdtd::geometry::Scene &sc = fdtd::geometry::scene();            ///< what the builders draw into
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::geometry; // transitional: module state, geometry/state.h

/**
 * @brief 
 * 
 * @param centerx 
 * @param centery_bulk 
 * @param centerz 
 * @param w 
 * @param t 
 * @param l1_r 
 * @param g 
 * @param l2_t 
 * @param index 
 * @param middle_layer 
 * @param antenna_shift 
 * @param l1_l 
 * @param l2_b 
 * @param alpha 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_lh_butterfly(const int &centerx, const int &centery_bulk, const int &centerz,
					   const double &w, const double &t, const double &/*l1_r*/, const double &g, const double &l2_t, const int &index, const double &middle_layer, const double &antenna_shift, const double &/*l1_l*/, const double &l2_b, const double &alpha)
{
	double beta_loc;
	double r;
	double r_arc;
	double centre_x_loc;
	double centre_z_loc;
	double beta_prev;
	r = double(corner_radius);
	r_arc = double(corner_radius);
	double cx_arc;
	double cz_arc;
	double cx_mono_prev;
	double cz_mono_prev;
	double cx_mono;
	double cz_mono;
	double delta = r * sin(PI / 8) / sin(3 * PI / 8);
	double delta_arc = r_arc * sin(PI / 8) / sin(3 * PI / 8);
	double coeff = 0.5 * sqrt(2) * g + delta; //0.5*sqrt(2)*g+r; //
	int i;

	double l_t;
	double l_r;
	double l_g; //l_t=l_b, l_r=l_l
	l_t = 0.5 * (l2_t + l2_t + w) - 0.5 * w + r - coeff; // length of the vertical monopoles
	l_r = 0.5 * (l2_b + l2_b + w) - 0.5 * w + r - coeff; // length of the horizontal monopoles
	l_g = sqrt(2) * (w - r + coeff) + 2 * r;			 //sqrt(2)*r; // length of the lateral monopoles that form the gap

	//dx
	beta_loc = double(0);
	centre_x_loc = 0.5 * w + coeff - r + 0.5 * l_r;
	centre_z_loc = 0;
	fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, l_r, t, w, index, middle_layer, antenna_shift, alpha, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

	//up
	beta_loc = double(0);
	centre_x_loc = 0;
	centre_z_loc = 0.5 * w + coeff - r + 0.5 * l_t;
	fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, l_t, index, middle_layer, antenna_shift, alpha, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

	//sx
	beta_loc = double(0);
	centre_x_loc = -(0.5 * w + coeff - r + 0.5 * l_r);
	centre_z_loc = 0;
	fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, l_r, t, w, index, middle_layer, antenna_shift, alpha, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

	//down
	beta_loc = double(0);
	centre_x_loc = 0;
	centre_z_loc = -(0.5 * w + coeff - r + 0.5 * l_t);
	fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, l_t, index, middle_layer, antenna_shift, alpha, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

	//g1
	beta_loc = double(-PI / 4);
	centre_x_loc = 0.25 * sqrt(2) * (w + g);
	centre_z_loc = 0.25 * sqrt(2) * (w + g);
	fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, l_g, t, w, index, middle_layer, antenna_shift, alpha, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

	//g2
	beta_loc = double(-PI / 4);
	centre_x_loc = -0.25 * sqrt(2) * (w + g);
	centre_z_loc = -0.25 * sqrt(2) * (w + g);
	fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, l_g, t, w, index, middle_layer, antenna_shift, alpha, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

	int angolo_tot = 45;
	int num_steps = 15;
	int step = angolo_tot / num_steps;

	//arc 1
	beta_prev = double(0);										  //valore iniziale
	cx_arc = (sqrt(2) - 0.5) * w + 0.5 * sqrt(2) * g + delta_arc; //centro arco
	cz_arc = 0.5 * w + r_arc;									  //centro arco
	cx_mono_prev = cx_arc;										  //centro_x monopolo iniziale
	cz_mono_prev = 0;											  //centro_z monopolo iniziale

	for (i = 0; i <= num_steps; i++)
	{
		beta_loc = beta_prev - (i * step * PI / 180);
		cx_mono = cx_arc + (cx_mono_prev - cx_arc) * cos(i * step * PI / 180) + (cz_mono_prev - cz_arc) * sin(i * step * PI / 180);
		cz_mono = cz_arc - (cx_mono_prev - cx_arc) * sin(i * step * PI / 180) + (cz_mono_prev - cz_arc) * cos(i * step * PI / 180);
		fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, w, index, middle_layer, antenna_shift, alpha, beta_loc, r, cx_mono, cz_mono, r, r);
	}

	//arc2
	beta_prev = double(0);										  //valore iniziale
	cz_arc = (sqrt(2) - 0.5) * w + 0.5 * sqrt(2) * g + delta_arc; //centro arco
	cx_arc = 0.5 * w + r_arc;									  //centro arco
	cx_mono_prev = 0;											  //centro_x monopolo iniziale
	cz_mono_prev = cz_arc;										  //centro_z monopolo iniziale

	for (i = 0; i <= num_steps; i++)
	{
		beta_loc = beta_prev - (i * step * PI / 180);
		cx_mono = cx_arc + (cx_mono_prev - cx_arc) * cos(-i * step * PI / 180) + (cz_mono_prev - cz_arc) * sin(-i * step * PI / 180);
		cz_mono = cz_arc - (cx_mono_prev - cx_arc) * sin(-i * step * PI / 180) + (cz_mono_prev - cz_arc) * cos(-i * step * PI / 180);
		fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, w, index, middle_layer, antenna_shift, alpha, beta_loc, r, cx_mono, cz_mono, r, r);
	}

	//arc 3
	beta_prev = double(0);											 //valore iniziale
	cx_arc = -((sqrt(2) - 0.5) * w + 0.5 * sqrt(2) * g + delta_arc); //centro arco
	cz_arc = -(0.5 * w + r_arc);									 //centro arco
	cx_mono_prev = cx_arc;											 //centro_x monopolo iniziale
	cz_mono_prev = 0;												 //centro_z monopolo iniziale

	for (i = 0; i <= num_steps; i++)
	{
		beta_loc = beta_prev - (i * step * PI / 180);
		cx_mono = cx_arc + (cx_mono_prev - cx_arc) * cos(i * step * PI / 180) + (cz_mono_prev - cz_arc) * sin(i * step * PI / 180);
		cz_mono = cz_arc - (cx_mono_prev - cx_arc) * sin(i * step * PI / 180) + (cz_mono_prev - cz_arc) * cos(i * step * PI / 180);
		fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, w, index, middle_layer, antenna_shift, alpha, beta_loc, r, cx_mono, cz_mono, r, r);
	}

	//arc4
	beta_prev = double(0);											 //valore iniziale
	cz_arc = -((sqrt(2) - 0.5) * w + 0.5 * sqrt(2) * g + delta_arc); //centro arco
	cx_arc = -(0.5 * w + r_arc);									 //centro arco
	cx_mono_prev = 0;												 //centro_x monopolo iniziale
	cz_mono_prev = cz_arc;											 //centro_z monopolo iniziale

	for (i = 0; i <= num_steps; i++)
	{
		beta_loc = beta_prev - (i * step * PI / 180);
		cx_mono = cx_arc + (cx_mono_prev - cx_arc) * cos(-i * step * PI / 180) + (cz_mono_prev - cz_arc) * sin(-i * step * PI / 180);
		cz_mono = cz_arc - (cx_mono_prev - cx_arc) * sin(-i * step * PI / 180) + (cz_mono_prev - cz_arc) * cos(-i * step * PI / 180);
		fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, w, index, middle_layer, antenna_shift, alpha, beta_loc, r, cx_mono, cz_mono, r, r);
	}

	return true;
}

/**
 * @brief Create a Round Tripole object
 * 
 * @param centerx 
 * @param centery_bulk 
 * @param centerz 
 * @param w 
 * @param t 
 * @param l1_r 
 * @param g 
 * @param l2_t 
 * @param index 
 * @param middle_layer 
 * @param antenna_shift 
 * @param l1_l 
 * @param l2_b 
 * @param alpha 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_round_tripole(const int &centerx, const int &centery_bulk, const int &centerz,
						const double &w, const double &t, const double &l1_r, const double &g, const double &l2_t, const int &index, const double &middle_layer, const double &antenna_shift, const double &l1_l, const double &l2_b, const double &alpha)
{
	double beta_loc;
	double r;
	double centre_x_loc;
	double centre_z_loc;
	r = double(corner_radius);

	//dx

	//up
	beta_loc = double(0);
	centre_x_loc = 0;
	centre_z_loc = (0.5 * l2_t + l1_l + 0.5 * l1_r);
	fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, l2_t, index, middle_layer, antenna_shift, alpha - PI / 4, beta_loc, r, centre_x_loc, centre_z_loc, r, r);
	//centre
	beta_loc = double(0);
	centre_x_loc = 0;
	centre_z_loc = 0;
	fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, l1_r, index, middle_layer, antenna_shift, alpha - PI / 4, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

	//sx

	//down
	beta_loc = double(0);
	centre_x_loc = 0;
	centre_z_loc = -(0.5 * l1_r + g + 0.5 * l2_b);
	fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, l2_b, index, middle_layer, antenna_shift, alpha - PI / 4, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

	////g1
	//
	////g2

	return true;
}

/**
 * @brief Create a Round Hybrid object
 * 
 * @param centerx 
 * @param centery_bulk 
 * @param centerz 
 * @param w 
 * @param t 
 * @param l1_r 
 * @param g 
 * @param l2_t 
 * @param index 
 * @param middle_layer 
 * @param antenna_shift 
 * @param l1_l 
 * @param l2_b 
 * @param alpha 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_round_hybrid(const int &centerx, const int &centery_bulk, const int &centerz,
					   const double &w, const double &t, const double &/*l1_r*/, const double &g, const double &l2_t, const int &index, const double &middle_layer, const double &antenna_shift, const double &/*l1_l*/, const double &l2_b, const double &alpha)
{
	int code;
	double beta_loc;
	double r;
	double centre_x_loc;
	double centre_z_loc;
	r = double(corner_radius);

	double l_t;						   //, l_r, l_g;  //l_t=l_b, l_r=l_l
	l_t = 0.5 * (l2_t + l2_b + w - g); // length of the vertical monopoles

	if (sc.draw)
	{
		//dx

		//centre
		if (g > 0)
		{
			beta_loc = double(0);
			centre_x_loc = 0;
			centre_z_loc = 0;
			code = sc.gap_material;
			fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, g, code, middle_layer, antenna_shift, alpha - PI / 4, beta_loc, double(0.E-9), centre_x_loc, centre_z_loc, double(0.E-9), double(0.E-9));
		}

		code = index;

		//up
		beta_loc = double(0);
		centre_x_loc = 0;
		centre_z_loc = 0.5 * (l_t + g);
		fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, l_t, code, middle_layer, antenna_shift, alpha - PI / 4, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

		//sx

		//down
		beta_loc = double(0);
		centre_x_loc = 0;
		centre_z_loc = -(0.5 * (l_t + g));
		fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, l_t, code, middle_layer, antenna_shift, alpha - PI / 4, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

		////g1
		//
		////g2
	}
	return true;
}

/**
 * @brief Create a Background object
 * 
 * @param centery_bulk 
 * @param centerx 
 * @param centerz 
 * @param middle_layer 
 * @param t 
 * @param antenna_shift 
 * @param index 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_background(const int &centery_bulk, const int &/*centerx*/, const int &/*centerz*/, const double &middle_layer, const double &/*t*/, const double &/*antenna_shift*/, const int &/*index*/)
{
	int i;
	int j;
	int k;

	int deep_y = centery_bulk - ceil(middle_layer / lat.d[1]);

	for (i = subdomain.output_start[0]; i < subdomain.output_start[0] + subdomain.output_size[0]; ++i)
	{
		for (j = subdomain.output_start[1]; j < subdomain.output_start[1] + subdomain.output_size[1]; ++j)
		{
			for (k = subdomain.output_start[2]; k < subdomain.output_start[2] + subdomain.output_size[2]; ++k)
			{
				if (j < deep_y) //<=
				{
					matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[0];
					matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[0];
					matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[0];

				}
				else if (j >= deep_y && j < centery_bulk) // >   <=
				{
					matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[1];
					matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[1];
					matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[1];

				}
				else if (j >= centery_bulk) // >
				{
					matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[2];
					matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[2];
					matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[2];
				}
			}
		}
	}

	if (fdtd::grid::world_rank == ROOTPROCESS)
	{
		for (j = 0; j <= lat.inner(1); ++j)
		{
			if (j < deep_y)
				background_y[j + j_offset_uni] = sc.medium[0];
			else if (j >= deep_y && j < centery_bulk)
				background_y[j + j_offset_uni] = sc.medium[1];
			else if (j >= centery_bulk)
				background_y[j + j_offset_uni] = sc.medium[2];
		}
	}

	return true;
}

/**
 * @brief 
 * 
 * @param centerx 
 * @param centery_bulk 
 * @param centerz 
 * @param SphereRadius 
 * @param ShellRadius 
 * @param index 
 * @param indexShell 
 * @param antenna_shift 
 * @return true 
 * @return false 
 */

/**
 * @brief 
 * 
 * @param centerx 
 * @param centery_bulk 
 * @param centerz 
 * @param w 
 * @param t 
 * @param l1_r 
 * @param g 
 * @param l2_t 
 * @param index 
 * @param middle_layer 
 * @param antenna_shift 
 * @param l1_l 
 * @param l2_b 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_l_antenna(const int &centerx, const int &centery_bulk, const int &centerz,
					 const double &w, const double &t, const double &l1_r, const double &/*g*/, const double &l2_t, const int &index, const double &/*middle_layer*/, const double &antenna_shift, const double &/*l1_l*/, const double &/*l2_b*/)
{
	int i;
	int j;
	int k;
	double x;
	double xf;
	double y;
	double yf;
	double z;
	double zf;
	double rEx;
	double rEy;
	double rEz;
	int centery = centery_bulk - int(antenna_shift / lat.d[1]); //- int(t/DY/2)

	for (i = subdomain.output_start[0]; i < subdomain.output_start[0] + subdomain.output_size[0]; ++i)

	{

		x = (i - centerx) * lat.d[0];
		xf = x + double(0.5) * lat.d[0];

		for (j = subdomain.output_start[1]; j < subdomain.output_start[1] + subdomain.output_size[1]; ++j)
		{

			y = (j - centery) * lat.d[1];
			yf = y + double(0.5) * lat.d[1];

			for (k = subdomain.output_start[2]; k < subdomain.output_start[2] + subdomain.output_size[2]; ++k)
			{

				z = (k - centerz) * lat.d[2];
				zf = z + double(0.5) * lat.d[2];

				if (sc.draw)
				{

					// vertical dipole
					x = (i - centerx) * lat.d[0] + double(0.5) * l1_r - double(0.5) * w;
					xf = x + double(0.5) * lat.d[0];

					if ((xf <= double(0.5) * w && xf >= -double(0.5) * w) && (y < 0 && y >= -t) && (z <= double(0.5) * l2_t && z >= -double(0.5) * l2_t))

						matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					if ((x <= double(0.5) * w && x >= -double(0.5) * w) && (yf <= 0 && yf > -t) && (z <= double(0.5) * l2_t && z >= -double(0.5) * l2_t))

						matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					if ((x <= double(0.5) * w && x >= -double(0.5) * w) && (y < 0 && y >= -t) && (zf <= double(0.5) * l2_t && zf >= -double(0.5) * l2_t))

						matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					//horizontal dipole
					x = (i - centerx) * lat.d[0];
					xf = x + double(0.5) * lat.d[0];
					z = (k - centerz) * lat.d[2] + double(0.5) * l2_t - double(0.5) * w;
					zf = z + double(0.5) * lat.d[2];

					if ((xf <= double(0.5) * l1_r && xf >= -double(0.5) * l1_r) && (y < 0 && y >= -t) && (z <= double(0.5) * w && z >= -double(0.5) * w))

						matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					if ((x <= double(0.5) * l1_r && x >= -double(0.5) * l1_r) && (yf <= 0 && yf > -t) && (z <= double(0.5) * w && z >= -double(0.5) * w))

						matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					if ((x <= double(0.5) * l1_r && x >= -double(0.5) * l1_r) && (y < 0 && y >= -t) && (zf <= double(0.5) * w && zf >= -double(0.5) * w))

						matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					// cilindro in alto

					x = (i - centerx) * lat.d[0] + double(0.5) * l1_r - double(0.5) * w;

					xf = x + double(0.5) * lat.d[0];

					z = (k - centerz) * lat.d[2] - double(0.5) * l2_t;

					zf = z + double(0.5) * lat.d[2];

					rEx = double(sqrt(double(xf * xf + z * z)));

					rEy = double(sqrt(double(x * x + z * z)));

					rEz = double(sqrt(double(x * x + zf * zf)));

					if ((rEx <= double(0.5) * w) && (y < 0 && y >= -t))
						matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;

					if ((rEy <= double(0.5) * w) && (yf <= 0 && yf > -t))
						matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;

					if ((rEz <= double(0.5) * w) && (y < 0 && y >= -t))
						matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;

					// cilindro a destra

					x = (i - centerx) * lat.d[0] - double(0.5) * l1_r;

					xf = x + double(0.5) * lat.d[0];

					z = (k - centerz) * lat.d[2] + double(0.5) * l2_t - double(0.5) * w;

					zf = z + double(0.5) * lat.d[2];

					rEx = double(sqrt(double(xf * xf + z * z)));

					rEy = double(sqrt(double(x * x + z * z)));

					rEz = double(sqrt(double(x * x + zf * zf)));

					if ((rEx <= double(0.5) * w) && (y < 0 && y >= -t))
						matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;

					if ((rEy <= double(0.5) * w) && (yf <= 0 && yf > -t))
						matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;

					if ((rEz <= double(0.5) * w) && (y < 0 && y >= -t))
						matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
				}
			}
		}
	}
	return true;
}

