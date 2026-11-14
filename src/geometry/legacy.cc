/**
 * @file geometry/legacy.cc
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
#include "core/constants.h"
#include "grid/subdomain.h"
#include "materials/material_map.h"
#include "materials/codes.h"
#include "grid/state.h"
#include "materials/state.h"
#include "io/text.h"
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
using fdtd::materials::nl_code;

namespace
{
const fdtd::geometry::Scene &sc = fdtd::geometry::scene();            ///< what the builders draw into
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace

/**
 * @brief Create a Xshaped object
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
 * @param r 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_x_shaped(const int &centerx, const int &centery_bulk, const int &centerz,
				   const double &w, const double &t, const double &l1_r, const double &g, const double &l2_t, const int &index, const double &middle_layer, const double &antenna_shift, const double &/*l1_l*/, const double &/*l2_b*/, const double &alpha, const double &r)
{
	double beta_loc, centre_x_loc, centre_z_loc;

	double L1 = l2_t;
	double w1 = w; // fixed monopole
	double w2 = g; // rotating monopole
	double L2 = l1_r;

	//fixed monopole

	beta_loc = double(0);
	centre_x_loc = 0;
	centre_z_loc = 0;
	fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, L1, t, w1, index, middle_layer, antenna_shift, 0 + PI / 4, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

	//rotating monopole
	beta_loc = -0.5 * PI + alpha; // double(0);
	centre_x_loc = 0;
	centre_z_loc = 0;
	fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w2, t, L2, index, middle_layer, antenna_shift, 0 + PI / 4, beta_loc, r, centre_x_loc, centre_z_loc, r, r);

	////g1
	//
	////g2

	return true;
}

/**
 * @brief Create a Hole object
 * 
 * @param centerx 
 * @param centery_bulk 
 * @param centerz 
 * @param Radius 
 * @param t 
 * @param antenna_shift 
 * @param index 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_hole(const int &centerx, const int &centery_bulk, const int &centerz, const double &Radius, const double &t, const double &antenna_shift, const int &index)
{
	int i, j, k;
	double x, xf, y, yf, z, zf, rEx, rEy, rEz;
	int centery = centery_bulk - int(antenna_shift / lat.d[1]);

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
					rEx = double(sqrt(double(xf * xf + z * z)));
					rEy = double(sqrt(double(z * z + x * x)));
					rEz = double(sqrt(double(x * x + zf * zf)));

					if (sc.staircase)
					{
						if (rEx > Radius && rEy > Radius && rEz > Radius && (y <= 0 && y >= -t))
						{
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
						}
					}
					else
					{
						if (rEx > Radius && (y <= 0 && y >= -t))
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
						if (rEy > Radius && (yf <= 0 && yf >= -t))
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
						if (rEz > Radius && (y <= 0 && y >= -t))
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
					}
				}
			}
		}
	}

	return true;
}

/**
 * @brief Create a Lateral Dipole object
 * 
 * @param centerx 
 * @param centery_bulk 
 * @param centerz 
 * @param w 
 * @param t 
 * @param l_t 
 * @param g 
 * @param l_b 
 * @param index 
 * @param middle_layer 
 * @param antenna_shift 
 * @param l_g 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_lateral_dipole(const int &centerx, const int &centery_bulk, const int &centerz,
						 const double &w, const double &t, const double &l_t, const double &g, const double &l_b, const int &index, const double &middle_layer, const double &antenna_shift, const double &l_g)
{
	int i, j, k;
	double x, xf, y, yf, z, zf, rEx, rEy, rEz;
	int centery = centery_bulk - int(antenna_shift / lat.d[1]); //- int(t/DY/2)
	int deep_y = centery_bulk - ceil(middle_layer / lat.d[1]);

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

					// vertical monopole top
					if ((xf <= -double(0.5) * g && xf >= -w - double(0.5) * g) && (y < 0 && y >= -t) && (z <= l_t - l_g / 2 && z >= -l_g / 2))

						matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					if ((x <= -double(0.5) * g && x >= -w - double(0.5) * g) && (yf <= 0 && yf > -t) && (z <= l_t - l_g / 2 && z >= -l_g / 2))

						matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					if ((x <= -double(0.5) * g && x >= -w - double(0.5) * g) && (y < 0 && y >= -t) && (zf <= l_t - l_g / 2 && zf >= -l_g / 2))

						matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					// vertical monopole bottom
					if ((xf <= w + double(0.5) * g && xf >= double(0.5) * g) && (y < 0 && y >= -t) && (z <= l_g / 2 && z >= -l_b + l_g / 2))

						matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					if ((x <= w + double(0.5) * g && x >= double(0.5) * g) && (yf <= 0 && yf > -t) && (z <= l_g / 2 && z >= -l_b + l_g / 2))

						matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					if ((x <= w + double(0.5) * g && x >= double(0.5) * g) && (y < 0 && y >= -t) && (zf <= l_g / 2 && zf >= -l_b + l_g / 2))

						matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					//
					//// rettangolo laterale a 45 gradi

					// gap dx

					if ((xf <= g + w && xf >= 0) && (y < 0 && y >= -t) && z > (-xf + (g + l_g) / 2))
					{
						if (j < deep_y)
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[0];
						else if (j >= deep_y && j < centery_bulk)
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[1];
						else if (j >= centery_bulk)
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[2];
					}

					if ((x <= g + w && x >= 0) && (yf < 0 && yf >= -t) && z > (-x + (g + l_g) / 2))
					{
						if (j < deep_y)
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[0];
						else if (j >= deep_y && j < centery_bulk)
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[1];
						else if (j >= centery_bulk)
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[2];
					}

					if ((x <= g + w && x >= 0) && (y < 0 && y >= -t) && zf > (-x + (g + l_g) / 2))
					{
						if (j < deep_y)
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[0];
						else if (j >= deep_y && j < centery_bulk)
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[1];
						else if (j >= centery_bulk)
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[2];
					}

					// gap sx

					if ((xf <= 0 && xf >= -g - w) && (y < 0 && y >= -t) && z < (-xf - (g + l_g) / 2))
					{
						if (j < deep_y)
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[0];
						else if (j >= deep_y && j < centery_bulk)
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[1];
						else if (j >= centery_bulk)
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[2];
					}

					if ((x <= 0 && x >= -g - w) && (yf < 0 && yf >= -t) && z < (-x - (g + l_g) / 2))
					{
						if (j < deep_y)
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[0];
						else if (j >= deep_y && j < centery_bulk)
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[1];
						else if (j >= centery_bulk)
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[2];
					}

					if ((x <= 0 && x >= -g - w) && (y < 0 && y >= -t) && zf < (-x - (g + l_g) / 2))
					{
						if (j < deep_y)
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[0];
						else if (j >= deep_y && j < centery_bulk)
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[1];
						else if (j >= centery_bulk)
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = sc.medium[2];
					}

					// cilindro in alto

					x = (i - centerx) * lat.d[0] + (w + g) / 2;

					xf = x + double(0.5) * lat.d[0];

					z = (k - centerz) * lat.d[2] - (l_t - l_g / 2);

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

					//cilindro in basso

					x = (i - centerx) * lat.d[0] - (w + g) / 2;

					xf = x + double(0.5) * lat.d[0];

					z = (k - centerz) * lat.d[2] + (l_b - l_g / 2);

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

					x = (i - centerx) * lat.d[0];

					xf = x + double(0.5) * lat.d[0];

					z = (k - centerz) * lat.d[2];

					zf = z + double(0.5) * lat.d[2];
				}
			}
		}
	}
	return true;
}

/**
 * @brief Create a Cube N L object
 * 
 * @param centery_bulk 
 * @param middle_layer 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_nonlinear_cube(const int &centery_bulk, const double &middle_layer)
{
	int i, j, k;
	int code;
	double x, xf, y, yf, z, zf, x_lim, y_lim, z_lim; // y_lim_min, y_lim_max;
	int deep_y = centery_bulk - ceil(middle_layer / lat.d[1]);
	x_lim = lat.d[0] * int((lat.cells[0] - 2 * lat.cpml[0] - 4 * lat.tfsf_inset[0]) / 2); //DX*lat.cells[0]; //
	y_lim = lat.d[1] * int((lat.cells[1] - 2 * lat.cpml[1] - 4 * lat.tfsf_inset[1]) / 2); //DY*lat.cells[1]; //
	z_lim = lat.d[2] * int((lat.cells[2] - 2 * lat.cpml[2] - 4 * lat.tfsf_inset[2]) / 2); //DZ*lat.cells[2]; //
	int centerx = int(lat.inner(0) / 2);
	int centery = int(lat.inner(1) / 2);
	int centerz = int(lat.inner(2) / 2);
	if (fdtd::grid::world_rank == ROOTPROCESS)
	{
		fdtd::io::log("geometry", "center NL cube: %d %d %d\n", centerx, centery, centerz);
		fdtd::io::log("geometry", "limits NL cube [nm]: %f %f %f\n", x_lim / double(1.E-9), y_lim / double(1.E-9), z_lim / double(1.E-9));
	}
	for (i = subdomain.output_start[0]; i < subdomain.output_start[0] + subdomain.output_size[0]; ++i)

	{

		x = (i - centerx) * lat.d[0];
		xf = x + double(0.5) * lat.d[0];

		for (j = subdomain.output_start[1]; j < subdomain.output_start[1] + subdomain.output_size[1]; ++j)
		{

			y = (j - centery) * lat.d[1];
			yf = y + double(0.5) * lat.d[1];

			if (j < deep_y)
				code = sc.medium[0] + nl_code;
			else if (j >= deep_y && j < centery_bulk)
				code = sc.medium[1] + nl_code;
			else if (j >= centery_bulk)
				code = sc.medium[2] + nl_code;

			for (k = subdomain.output_start[2]; k < subdomain.output_start[2] + subdomain.output_size[2]; ++k)
			{

				z = (k - centerz) * lat.d[2];
				zf = z + double(0.5) * lat.d[2];

				if ((xf <= x_lim && xf >= -x_lim) && (y <= y_lim && y >= -y_lim) && (z <= z_lim && z >= -z_lim))

					matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = code;

				if ((x <= x_lim && x >= -x_lim) && (yf <= y_lim && yf >= -y_lim) && (z <= z_lim && z >= -z_lim))

					matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = code;

				if ((x <= x_lim && x >= -x_lim) && (y <= y_lim && y >= -y_lim) && (zf <= z_lim && zf >= -z_lim))

					matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = code;
			}
		}
	}
	return true;
}

/**
 * @brief Create a Dipole Round object
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
bool fdtd::geometry::create_dipole_round(const int &centerx, const int &centery_bulk, const int &centerz,
					   const double &w, const double &t, const double &l1_r, const double &g, const double &/*l2_t*/, const int &index, const double &/*middle_layer*/, const double &antenna_shift, const double &l1_l, const double &/*l2_b*/, const double &alpha)
{
	int i, j, k;
	double x, xf, y, yf, z, zf, rEx, rEy, rEz, _x, _xf, _z, _zf;
	int centery = centery_bulk - int(antenna_shift / lat.d[1]); //- int(t/DY/2)
	double theta_ref = double(0);
	double beta = alpha - theta_ref;

	for (i = subdomain.output_start[0]; i < subdomain.output_start[0] + subdomain.output_size[0]; ++i)

	{

		_x = (i - centerx) * lat.d[0];
		_xf = _x + double(0.5) * lat.d[0];

		for (j = subdomain.output_start[1]; j < subdomain.output_start[1] + subdomain.output_size[1]; ++j)
		{

			y = (j - centery) * lat.d[1];
			yf = y + double(0.5) * lat.d[1];

			for (k = subdomain.output_start[2]; k < subdomain.output_start[2] + subdomain.output_size[2]; ++k)
			{

				_z = (k - centerz) * lat.d[2];
				_zf = _z + double(0.5) * lat.d[2];

				x = _x * cos(beta) + _z * sin(beta);
				xf = _xf * cos(beta) + _z * sin(beta);
				z = -_x * sin(beta) + _z * cos(beta);
				zf = -_x * sin(beta) + _zf * cos(beta);

				if (sc.draw)
				{

					// dipole

					if (sc.staircase)
					{
						if (((xf <= l1_r && xf >= 0.5 * g) || (xf <= -0.5 * g && xf >= -l1_l)) && (yf < 0 && yf >= -t) && (zf <= double(0.5) * w && zf >= -double(0.5) * w))
						{
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);
						}
					}

					else
					{
						if (((xf <= l1_r && xf >= 0.5 * g) || (xf <= -0.5 * g && xf >= -l1_l)) && (y < 0 && y >= -t) && (z <= double(0.5) * w && z >= -double(0.5) * w))

							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

						if (((x <= l1_r && x >= 0.5 * g) || (x <= -0.5 * g && x >= -l1_l)) && (yf <= 0 && yf > -t) && (z <= double(0.5) * w && z >= -double(0.5) * w))

							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

						if (((xf <= l1_r && xf >= 0.5 * g) || (xf <= -0.5 * g && xf >= -l1_l)) && (y < 0 && y >= -t) && (zf <= double(0.5) * w && zf >= -double(0.5) * w))

							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);
					}
					// cilindro a destra

					x = _x * cos(beta) + _z * sin(beta) - l1_r;
					xf = _xf * cos(beta) + _z * sin(beta) - l1_r;
					z = -_x * sin(beta) + _z * cos(beta);
					zf = -_x * sin(beta) + _zf * cos(beta);

					rEx = double(sqrt(double(xf * xf + z * z)));

					rEy = double(sqrt(double(x * x + z * z)));

					rEz = double(sqrt(double(x * x + zf * zf)));

					if (sc.staircase)
					{
						if ((rEx <= double(0.5) * w) && (rEy <= double(0.5) * w) && (rEz <= double(0.5) * w) && (yf <= 0 && yf > -t))
						{
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
						}
					}
					else
					{
						if ((rEx <= double(0.5) * w) && (y < 0 && y >= -t))
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;

						if ((rEy <= double(0.5) * w) && (yf <= 0 && yf > -t))
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;

						if ((rEz <= double(0.5) * w) && (y < 0 && y >= -t))
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
					}

					//cilindro a sinistra

					x = _x * cos(beta) + _z * sin(beta) + l1_l;
					xf = _xf * cos(beta) + _z * sin(beta) + l1_l;
					z = -_x * sin(beta) + _z * cos(beta);
					zf = -_x * sin(beta) + _zf * cos(beta);

					rEx = double(sqrt(double(xf * xf + z * z)));

					rEy = double(sqrt(double(x * x + z * z)));

					rEz = double(sqrt(double(x * x + zf * zf)));

					if (sc.staircase)
					{
						if ((rEx <= double(0.5) * w) && (rEy <= double(0.5) * w) && (rEz <= double(0.5) * w) && (yf <= 0 && yf > -t))
						{
							matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
							matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
							matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = (short)index;
						}
					}

					else
					{
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
	}

	return true;
}
