/**
 * @file geometry/basic.cc
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
#include "sources/signal.h"
#include "grid/subdomain.h"
#include "grid/state.h"
#include "io/text.h"
#include "geometry/state.h"
#include "geometry/builders.h"
#include "geometry/render.h"
#include "grid/lattice.h"
#include "geometry/scene.h"
#include "geometry/block.h"
#include "geometry/shapes.h"
using fdtd::constants::PI;
using fdtd::grid::subdomain;
using fdtd::materials::matind_x;
using fdtd::materials::matind_y;
using fdtd::materials::matind_z;

namespace
{
const fdtd::geometry::Scene &sc = fdtd::geometry::scene();            ///< what the builders draw into
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::geometry; // transitional: module state, geometry/state.h

/**
 * @brief Rounded box, "rounded monopole": a b x h x t block with every edge and corner rounded
 *
 * @verbatim
      top view, local x-z,                      side view, local x-y
        z ^       r_c                             y ^
          |   .---------.   ---  h/2               0 |    .-----------.      <- edges rounded by r_f
          |  /           \                           |   /             \
          | |      c      |    h                     |  |       c       |    t
          | |   (cx,cz)   |                          |   \             /
          |  \           /                        -t |    `-----------'      <- edges rounded by r_n
          |   `---------'   --- -h/2                 +----------------------> x
          |   |<-- b -->|                                 |<--- b --->|
          +-------------------> x
   @endverbatim
 *
 * The block spans y in [-t, 0] below the reference height centery_bulk dy - antenna_shift
 * antenna_shift is not rounded to cells here. Its rounding:
 *   - the four vertical edges, along y: radius r_c (`r`);
 *   - the edges of the top face (y = 0): radius r_f (`round_far`), of the bottom face: r_n (`round_near`);
 *   - the eight corners: spheres of radius r_c.
 * Each slice at height y is drawn as boxes and discs, a sliced solid.
 *
 * Placement: the block is drawn about the local centre (centre_x_loc, centre_z_loc), turned by
 * beta_loc about that centre, then the whole drawing is turned by alpha - pi/4 about the centre cell:
 * @verbatim
        global = R(alpha - pi/4) [ R(beta_loc) block + (centre_x_loc, centre_z_loc) ]
   @endverbatim
 * Staircasing: all three Yee sample points inside.
 *
 * @param b, h, t      block size along local x, local z and y [m]
 * @param index        material
 * @param middle_layer unused
 * @param r            radius of the vertical edges and corner spheres (r_c) [m]
 */
bool fdtd::geometry::create_round_monopole(const int &centerx, const int &centery_bulk, const int &centerz,
						 const double &b, const double &t, const double &h, const int &index, const double &middle_layer, const double &antenna_shift, const double &alpha, const double &beta_loc, const double &r, const double &centre_x_loc, const double &centre_z_loc, const double &round_near, const double &round_far)
{
	(void)middle_layer;
	if (!sc.draw)
		return true;
	using namespace geom;
	const double r_c = r;
	const double r_f = round_far;
	const double T = t;
	const double r_n = double(round_near);
	const double ex = 0.5 * b - r_c;
	const double ez = 0.5 * h - r_c; // centres of the vertical edges

	Solid block = extrude(rounded_rect(0, 0, b, h, r_c), -T + r_n, -r_f) + extrude(box(-ex, -ez, ex, ez), -T, 0);
	for (auto [rad, yc] : {std::pair{r_f, -r_f}, std::pair{r_n, -T + r_n}}) // rounded edges of the top and bottom faces
		for (double s : {-1.0, 1.0})
			block = block + cylinder_z(s * (0.5 * b - rad), yc, rad, -ez, ez) + cylinder_x(yc, s * (0.5 * h - rad), rad, -ex, ex);
	for (double yc : {-r_c, -T + r_c}) // corner balls
		for (auto [sx, sz] : {std::pair{-1, -1}, {-1, 1}, {1, -1}, {1, 1}})
			block = block + geom::sphere(sx * ex, yc, sz * ez, r_c);

	const Solid placed = block.rotated(beta_loc, 0, 0).translated(centre_x_loc, -antenna_shift, centre_z_loc).rotated(alpha - 0.25 * PI, 0, 0);
	render_parts({{placed, index}}, centerx, centery_bulk, centerz);
	return true;
}

/**
 * @brief Two rounded blocks, CreateRoundMonopole, on the z axis, separated by the gap g
 * @image html round_dipole.png "[round_dipole] keys on the mesh, docs/geometry/Shapes.md"
 * Each block is w x l x t with radii round_centre, vertical edges, corners, round_far, top-face
 * edges, and round_near, bottom-face edges; see CreateRoundMonopole. l1_r, l1_l and middle_layer
 * are unused.
 */
bool fdtd::geometry::create_round_dipole(const int &centerx, const int &centery_bulk, const int &centerz,
					   const double &w, const double &t, const double &l1_r, const double &g, const double &l2_t, const int &index, const double &middle_layer, const double &antenna_shift, const double &l1_l, const double &l2_b, const double &alpha, const double &round_centre, const double &round_near, const double &round_far)
{
	(void)l1_r, (void)l1_l;
	for (auto [l, z] : {std::pair{l2_t, double(0.5) * (l2_t + g)}, std::pair{l2_b, -double(0.5) * (l2_b + g)}})
		fdtd::geometry::create_round_monopole(centerx, centery_bulk, centerz, w, t, l, index, middle_layer, antenna_shift, alpha - double(0.25 * PI), 0, round_centre, 0, z, round_near, round_far);
	return true;
}

/**
 * @brief Sphere of radius `SphereRadius` centred on a cell
 *
 * @image html sphere.png "[sphere] keys on the mesh, docs/geometry/Shapes.md"
 *
 * Staircasing: a cell is filled when all three Yee sample points are inside.
 *
 * @param centerx, centery_bulk, centerz  centre cell, inner coordinates
 * @param SphereRadius                    radius [m]
 * @param index                           material
 * @param antenna_shift                   shift of the centre towards -y [m], rounded down to whole cells
 */
bool fdtd::geometry::create_sphere(const int &centerx, const int &centery_bulk, const int &centerz,
				  const double &SphereRadius, const int &index, const double &antenna_shift)
{
	if (!sc.draw)
		return true;
	render_parts({{geom::sphere(0, 0, 0, SphereRadius), index}}, centerx, centery_bulk - int(antenna_shift / lat.d[1]), centerz);
	return true;
}

/**
 * @brief Create a Sphere 2 object
 * 
 * @param centerx 
 * @param centery_bulk 
 * @param centerz 
 * @param SphereRadius 
 * @param index 
 * @param antenna_shift 
 * @return true 
 * @return false 
 */
//overloaded to accept in input double coordinates instead of integer coordinates
/**
 * @brief Sphere whose centre is given in metres (x, z) and need not sit on a cell
 *
 * @image html sphere_dimer.png "[sphere_dimer] keys on the mesh, docs/geometry/Shapes.md"
 * Staircasing: all three Yee sample points inside.
 */
bool fdtd::geometry::create_sphere_2(const double &centerx, const int &centery_bulk, const double &centerz,
				   const double &SphereRadius, const int &index, const double &antenna_shift)
{
	if (!sc.draw)
		return true;
	render_parts({{geom::sphere(centerx, -antenna_shift, centerz, SphereRadius), index}}, 0, centery_bulk, 0);
	return true;
}

/**
 * @brief Rectangular monopole with its four vertical edges rounded
 * @image html rounded_monopole.png "[rounded_monopole] keys on the mesh, docs/geometry/Shapes.md"
 * The outline is a rounded rectangle: two rectangles plus four discs of radius r.
 * Staircasing: the cell centre must be inside.
 *
 * @param w, l_t  size along x and z [m]
 * @param t       thickness along y [m]
 * @param radius  corner radius r [m], 0 gives sharp corners
 */
bool fdtd::geometry::create_monopole_rounded_lateral(const int &centerx, const int &centery_bulk, const int &centerz, const double &w, const double &t, const double &l_t, const double &radius, const int &index, const double &antenna_shift)
{
	if (!sc.draw)
		return true;
	render_parts({{geom::extrude(geom::rounded_rect(0, 0, w, l_t, radius), -0.5 * t, 0.5 * t), index}},
				centerx, centery_bulk - int(antenna_shift / lat.d[1]), centerz, Staircase::CellCentre);
	return true;
}

/**
 * @brief Gap dipole: two rectangular arms along z, separated by a gap, of width w and thickness t
 *
 * @image html dipole.png "[dipole] keys on the mesh, docs/geometry/Shapes.md"
 *
 * Staircasing: the cell centre must be inside.
 *
 * @param centerx, centery_bulk, centerz  centre of the gap
 * @param w    arm width (x) [m]
 * @param t    thickness (y) [m]
 * @param l_t  end of the top arm, measured from the centre (z) [m]
 * @param g    gap [m]
 * @param l_b  end of the bottom arm, measured from the centre (-z) [m]
 * @param index          material
 * @param antenna_shift  shift towards -y [m], whole cells
 */
bool fdtd::geometry::create_dipole(const int &centerx, const int &centery_bulk, const int &centerz,
				  const double &w, const double &t, const double &l_t, const double &g, const double &l_b, const int &index, const double &antenna_shift)
{
	if (!sc.draw)
		return true;
	const geom::Shape arms = geom::box(-0.5 * w, 0.5 * g, 0.5 * w, l_t) + geom::box(-0.5 * w, -l_b, 0.5 * w, -0.5 * g);
	render_parts({{geom::extrude(arms, -0.5 * t, 0.5 * t), index}}, centerx, centery_bulk - int(antenna_shift / lat.d[1]), centerz,
				Staircase::CellCentre);
	return true;
}

/**
 * @brief Bowtie: two triangles along z, tips facing each other across the gap
 *
 * @image html bowtie.png "[bowtie] keys on the mesh, docs/geometry/Shapes.md"
 *
 * Thickness t along y, centred on the centre cell. Staircasing: the cell centre must be inside.
 *
 * @param w  width of the triangle bases [m]
 * @param t  thickness (y) [m]
 * @param l  total length, tip to tip of the outer bases (z) [m]
 * @param g  gap [m]
 */
bool fdtd::geometry::create_bowtie(const int &centerx, const int &centery_bulk, const int &centerz,
				  const double &w, const double &t, const double &l, const double &g, const int &index, const double &antenna_shift)
{
	if (!sc.draw)
		return true;
	const double m = l / w; // |x| <= |z| / m
	auto wing = [&](double sign) {
		const double z0 = sign * 0.5 * g;
		const double z1 = sign * 0.5 * l;
		return geom::polygon({{-0.5 * g / m, z0}, {0.5 * g / m, z0}, {0.5 * l / m, z1}, {-0.5 * l / m, z1}});
	};
	render_parts({{geom::extrude(wing(+1) + wing(-1), -0.5 * t, 0.5 * t), index}}, centerx, centery_bulk - int(antenna_shift / lat.d[1]), centerz,
				Staircase::CellCentre);
	return true;
}

/**
 * @brief Disc in the x-z plane, cylinder along y, of radius `Radius` and height `Height`
 *
 * @image html cylinder_3d.png "[cylinder_3d] keys on the mesh, docs/geometry/Shapes.md"
 *
 * The height is centred on the centre cell. Staircasing: all three Yee points inside.
 *
 * @param centerx, centery_bulk, centerz  centre cell
 * @param Radius                          [m]
 * @param Height                          [m], along y
 * @param antenna_shift                   shift towards -y [m], whole cells
 * @param index                           material
 */
bool fdtd::geometry::create_cylinder_xz(const int &centerx, const int &centery_bulk, const int &centerz, const double &Radius, const double &Height, const double &antenna_shift, const int &index)
{
	if (!sc.draw)
		return true;
	render_parts({{geom::extrude(geom::circle(0, 0, Radius), -0.5 * Height, 0.5 * Height), index}},
				centerx, centery_bulk - int(antenna_shift / lat.d[1]), centerz);
	return true;
}

/**
 * @brief Triangle in the x-z plane, base at -l/2, apex at +l/2
 *
 * @image html triangle.png "[triangle] keys on the mesh, docs/geometry/Shapes.md"
 *
 * Thickness t along y, centred on the centre cell. Staircasing: the cell centre must be inside
 * the original routine filled nothing when staircasing.
 *
 * @param w  base width [m]
 * @param t  thickness (y) [m]
 * @param l  height of the triangle (z) [m]
 * @param g  unused
 */
bool fdtd::geometry::create_triangle(const int &centerx, const int &centery_bulk, const int &centerz, const double &w, const double &t, const double &l, const double &g, const int &index, const double &antenna_shift)
{
	(void)g;
	if (!sc.draw)
		return true;
	const geom::Shape tri = geom::polygon({{-0.5 * w, -0.5 * l}, {0.5 * w, -0.5 * l}, {0, 0.5 * l}});
	render_parts({{geom::extrude(tri, -0.5 * t, 0.5 * t), index}}, centerx, centery_bulk - int(antenna_shift / lat.d[1]), centerz,
				Staircase::CellCentre);
	return true;
}

/**
 * @brief Create a Roatated Monopole object
 * 
 * @param centerx 
 * @param centery_bulk 
 * @param centerz 
 * @param w 
 * @param t 
 * @param l_t 
 * @param index 
 * @param middle_layer 
 * @param antenna_shift 
 * @param theta 
 * @param radius 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_rotated_monopole(const int &centerx, const int &centery_bulk, const int &centerz, const double &w, const double &t, const double &l_t, const int &index, const double &/*middle_layer*/, const double &antenna_shift, const double &theta, const double &radius)
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
	int centery = centery_bulk - ceil(antenna_shift / lat.d[1]); //- int(t/DY/2)

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
					if (((xf * cos(theta) - z * sin(theta)) <= double(0.5) * w && (xf * cos(theta) - z * sin(theta)) >= -double(0.5) * w) && (y <= t / 2.0 && y >= -t / 2.0) && ((z * cos(theta) + xf * sin(theta)) <= l_t / 2 && (z * cos(theta) + xf * sin(theta)) >= -l_t / 2))

						matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					if (((x * cos(theta) - z * sin(theta)) <= double(0.5) * w && (x * cos(theta) - z * sin(theta)) >= -double(0.5) * w) && (yf <= t / 2.0 && yf > -t / 2.0) && ((z * cos(theta) + x * sin(theta)) <= l_t / 2 && (z * cos(theta) + x * sin(theta)) >= -l_t / 2))

						matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);

					if (((x * cos(theta) - zf * sin(theta)) <= double(0.5) * w && (xf * cos(theta) - z * sin(theta)) >= -double(0.5) * w) && (y <= t / 2.0 && y >= -t / 2.0) && ((zf * cos(theta) + x * sin(theta)) <= l_t / 2 && (zf * cos(theta) + x * sin(theta)) >= -l_t / 2))

						matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = short(index);
				}
			}
		}
	}
	double dipolelength = 2.5 * radius;

	if ((sc.shape == shape::v_shape))
	{
		fdtd::geometry::create_round_monopole(centerx + (0.5 * w / lat.d[0]) * cos(theta) + (0.5 * l_t / lat.d[2]) * sin(theta), centery, centerz + (0.5 * l_t / lat.d[2]) * cos(theta) - (0.5 * w / lat.d[0]) * sin(theta), radius * sqrt(2), t, dipolelength, sc.medium[1], 0, -t / 2, theta + PI / 2, 0, 0, 0, 0, 0, 0);
		if (!((sc.shape == shape::v_shape) && theta > PI / 2.0))
			fdtd::geometry::create_round_monopole(centerx - (0.5 * w / lat.d[0]) * cos(theta) + (0.5 * l_t / lat.d[2]) * sin(theta), centery, centerz + (0.5 * l_t / lat.d[2]) * cos(theta) + (0.5 * w / lat.d[0]) * sin(theta), radius * sqrt(2), t, dipolelength, sc.medium[1], 0, -t / 2, -theta, 0, 0, 0, 0, 0, 0);
		fdtd::geometry::create_round_monopole(centerx + (0.5 * w / lat.d[0]) * cos(theta) - (0.5 * l_t / lat.d[2]) * sin(theta), centery, centerz - (0.5 * l_t / lat.d[2]) * cos(theta) - (0.5 * w / lat.d[0]) * sin(theta), radius * sqrt(2), t, dipolelength, sc.medium[1], 0, -t / 2, -theta, 0, 0, 0, 0, 0, 0);
		if (!((sc.shape == shape::v_shape) && theta < PI / 2.0))
			fdtd::geometry::create_round_monopole(centerx - (0.5 * w / lat.d[0]) * cos(theta) - (0.5 * l_t / lat.d[2]) * sin(theta), centery, centerz - (0.5 * l_t / lat.d[2]) * cos(theta) + (0.5 * w / lat.d[0]) * sin(theta), radius * sqrt(2), t, dipolelength, sc.medium[1], 0, -t / 2, theta + PI / 2, 0, 0, 0, 0, 0, 0);

		fdtd::geometry::create_cylinder_xz(centerx + (0.5 * w / lat.d[0] - radius / lat.d[0]) * cos(theta) + (0.5 * l_t / lat.d[2] - radius / lat.d[2]) * sin(theta), centery, centerz + (0.5 * l_t / lat.d[2] - radius / lat.d[2]) * cos(theta) - (0.5 * w / lat.d[0] - radius / lat.d[0]) * sin(theta), radius, t, 0, index);
		if (!((sc.shape == shape::v_shape) && theta > PI / 2.0))
			fdtd::geometry::create_cylinder_xz(centerx - (0.5 * w / lat.d[0] - radius / lat.d[0]) * cos(theta) + (0.5 * l_t / lat.d[2] - radius / lat.d[2]) * sin(theta), centery, centerz + (0.5 * l_t / lat.d[2] - radius / lat.d[2]) * cos(theta) + (0.5 * w / lat.d[0] - radius / lat.d[0]) * sin(theta), radius, t, 0, index);
		fdtd::geometry::create_cylinder_xz(centerx + (0.5 * w / lat.d[0] - radius / lat.d[0]) * cos(theta) - (0.5 * l_t / lat.d[2] - radius / lat.d[2]) * sin(theta), centery, centerz - (0.5 * l_t / lat.d[2] - radius / lat.d[2]) * cos(theta) - (0.5 * w / lat.d[0] - radius / lat.d[0]) * sin(theta), radius, t, 0, index);
		if (!((sc.shape == shape::v_shape) && theta < PI / 2.0))
			fdtd::geometry::create_cylinder_xz(centerx - (0.5 * w / lat.d[0] - radius / lat.d[0]) * cos(theta) - (0.5 * l_t / lat.d[2] - radius / lat.d[2]) * sin(theta), centery, centerz - (0.5 * l_t / lat.d[2] - radius / lat.d[2]) * cos(theta) + (0.5 * w / lat.d[0] - radius / lat.d[0]) * sin(theta), radius, t, 0, index);
	}

	if ((sc.shape == shape::l_shape))
	{
		fdtd::geometry::create_round_monopole(centerx + int(0.5 * w / lat.d[0]), centery, centerz + int(0.5 * l_t / lat.d[2]), radius *sqrt(2), t, dipolelength, sc.medium[1], 0, -t / 2, theta + PI / 2, 0, 0, 0, 0, 0, 0);
		if (!(centerx == sc.centre[0]))
			fdtd::geometry::create_round_monopole(centerx - int(0.5 * w / lat.d[0]), centery, centerz + int(0.5 * l_t / lat.d[2]), radius *sqrt(2), t, dipolelength, sc.medium[1], 0, -t / 2, -theta, 0, 0, 0, 0, 0, 0);
		if (!(centerz == sc.centre[2]))
			fdtd::geometry::create_round_monopole(centerx + int(0.5 * w / lat.d[0]), centery, centerz - int(0.5 * l_t / lat.d[2]), radius *sqrt(2), t, dipolelength, sc.medium[1], 0, -t / 2, -theta, 0, 0, 0, 0, 0, 0);
		fdtd::geometry::create_round_monopole(centerx - int(0.5 * w / lat.d[0]), centery, centerz - int(0.5 * l_t / lat.d[2]), radius *sqrt(2), t, dipolelength, sc.medium[1], 0, -t / 2, theta + PI / 2, 0, 0, 0, 0, 0, 0);

		fdtd::geometry::create_cylinder_xz(centerx + int(0.5 * w / lat.d[0] - radius / lat.d[0]), centery, centerz + int(0.5 * l_t / lat.d[2] - radius / lat.d[2]), radius, t, 0, index);
		if (!(centerx == sc.centre[0]))
			fdtd::geometry::create_cylinder_xz(centerx - int(0.5 * w / lat.d[0] - radius / lat.d[0]), centery, centerz + int(0.5 * l_t / lat.d[2] - radius / lat.d[2]), radius, t, 0, index);
		if (!(centerz == sc.centre[2]))
			fdtd::geometry::create_cylinder_xz(centerx + int(0.5 * w / lat.d[0] - radius / lat.d[0]), centery, centerz - int(0.5 * l_t / lat.d[2] - radius / lat.d[2]), radius, t, 0, index);
		fdtd::geometry::create_cylinder_xz(centerx - int(0.5 * w / lat.d[0] - radius / lat.d[0]), centery, centerz - int(0.5 * l_t / lat.d[2] - radius / lat.d[2]), radius, t, 0, index);
	}

	return true;
}

/**
 * @brief Isosceles triangle with all three corners rounded by r, overall size w x l, turned by theta
 * @image html rounded_triangle.png "[rounded_triangle] keys on the mesh, docs/geometry/Shapes.md"
 * w and l are the outer size of the rounded triangle, its bounding box, centred on the centre cell.
 * The drawing is turned by -theta, clockwise, about the centre cell. Thickness t along y, centred on
 * the reference height. Staircasing: the cell centre must be inside.
 *
 * Replaces an iterative search, 1000 x 1000 grid, cell tolerance, for the sharp triangle whose
 * rounded version has size w x l: the inner triangle gives it exactly.
 *
 * @param w, l    base width and height of the rounded triangle [m]
 * @param t       thickness (y) [m]
 * @param radius  corner radius r [m], at most min(w, l) / 2
 * @param theta   rotation [rad]
 */
bool fdtd::geometry::create_rounded_triangle(const int &centerx, const int &centery_bulk, const int &centerz, const double &w, const double &t, const double &l, const int &index, const double &antenna_shift, const double &radius, const double &theta)
{
	if (!sc.draw)
		return true;
	const double r = radius;
	const double hw = 0.5 * w - r;
	const double hl = 0.5 * l - r; // inner triangle
	const geom::Shape tri = geom::inflated_polygon({{-hw, -hl}, {hw, -hl}, {0, hl}}, r);
	render_parts({{geom::extrude(tri.rotated(-theta, 0, 0), -0.5 * t, 0.5 * t), index}},
				centerx, centery_bulk - int(antenna_shift / lat.d[1]), centerz, Staircase::CellCentre);
	return true;
}

/// area of CreateRoundedTriangle's drawing: inner area + inner perimeter * r + pi r^2
double fdtd::geometry::rounded_triangle_area(double w, double l, double r)
{
	const double b = w - 2 * r;
	const double h = l - 2 * r;
	return double(0.5) * b * h + (b + 2 * std::hypot(double(0.5) * b, h)) * r + PI * r * r;
}

/**
 * @brief Create a Layer object
 * 
 * @param start_y 
 * @param end_y 
 * @param medium 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_layer(const int &start_y, const int &end_y, const int &medium)
{
	// isotropic special case of CreateAnisLayer
	return fdtd::geometry::create_anisotropic_layer(start_y, end_y, medium, medium, medium);
}

/**
 * @brief Create a Anis Layer object
 * 
 * @param start_y 
 * @param end_y 
 * @param medium_x 
 * @param medium_y 
 * @param medium_z 
 * @return true 
 * @return false 
 */
bool fdtd::geometry::create_anisotropic_layer(const int &start_y, const int &end_y, const int &medium_x, const int &medium_y, const int &medium_z)
{
	int i;
	int j;
	int k;

	for (i = subdomain.output_start[0]; i < subdomain.output_start[0] + subdomain.output_size[0]; ++i)
	{
		for (j = subdomain.output_start[1]; j < subdomain.output_start[1] + subdomain.output_size[1]; ++j)
		{
			for (k = subdomain.output_start[2]; k < subdomain.output_start[2] + subdomain.output_size[2]; ++k)
			{
				if (j >= start_y && j < end_y)
				{
					matind_x[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = medium_x;
					matind_y[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = medium_y;
					matind_z[i - subdomain.output_start[0]][j - subdomain.output_start[1]][k - subdomain.output_start[2]] = medium_z;
				}
			}
		}
	}
	return true;
}

/**
 * @brief Infinite cylinder along z of radius `Radius`
 *
 * @image html cylinder.png "[cylinder] keys on the mesh, docs/geometry/Shapes.md"
 *
 * Staircasing: all three Yee points inside.
 *
 * @param centerx, centery  axis cell (x, y); centerz only sets the reference of z samples
 * @param Radius            [m]
 * @param index             material
 */
bool fdtd::geometry::create_cylinder(const int &centerx, const int &centery, const int &centerz, const double &Radius, const int &index, double length)
{
	if (!sc.draw)
		return true;
	// length 0: through the whole domain, the default span of geom::cylinder_z
	render_parts({{length > 0 ? geom::cylinder_z(0, 0, Radius, -0.5 * length, 0.5 * length) : geom::cylinder_z(0, 0, Radius), index}}, centerx, centery, centerz);
	return true;
}

/**
 * @brief Two coaxial cylinders along y, the propagation axis: the shell (radius_shell, shell_length) first, the core (radius_core, core_length) over it
 * @image html cylinder_shell.png "[cylinder_shell] keys on the mesh, docs/geometry/Shapes.md"
 * The two lengths are independent: a shell longer than the core closes the ends, a shorter one leaves the
 * ends of the core bare. A length of 0 runs that cylinder through the whole domain.
 * Staircasing: all three Yee points inside.
 *
 * @param radius_core, core_length    core [m]
 * @param radius_shell, shell_length  shell [m], outer radius
 * @param index, shell_index          materials of the core and of the shell
 */
bool fdtd::geometry::create_cylinder_shell(int centerx, int centery, int centerz, double radius_core, double core_length, double radius_shell, double shell_length, int index, int shell_index)
{
	if (!sc.draw)
		return true;
	// a disc in the x-z plane extruded along y; length 0 spans the domain, as cylinder_z
	auto cylinder = [](double radius, double length) {
		const double half = length > 0 ? 0.5 * length : 1e3;
		return radius > 0 ? geom::extrude(geom::circle(0, 0, radius), -half, half) : geom::Solid();
	};
	// in this order: the core overwrites the shell where they overlap
	render_parts({{cylinder(radius_shell, shell_length), shell_index}, {cylinder(radius_core, core_length), index}}, centerx, centery, centerz);
	return true;
}

/**
 * @brief The [block.*] boxes, in file order, about the structure centre; each block's copies are
 *        offset by multiples of its pitch. Staircasing: all three Yee points inside.
 */
void fdtd::geometry::create_blocks(const std::vector<Block>& blocks)
{
	std::vector<geom::Part> parts;
	for (const auto& b : blocks)
		for (int n = 0; n < b.count; ++n)
		{
			const double sx = n * b.pitch_x;
			const double sy = n * b.pitch_y;
			const double sz = n * b.pitch_z;
			parts.push_back({geom::extrude(geom::box(b.x_min + sx, b.z_min + sz, b.x_max + sx, b.z_max + sz), b.y_min + sy, b.y_max + sy), b.material});
		}
	if (!parts.empty())
		render_parts(parts, sc.centre[0], sc.centre[1], sc.centre[2]);
}

/**
 * @brief Elliptic disc in the x-z plane: semi-axes radius_x and radius_z, turned by theta, thickness along y
 * @image html ellipse.png "[ellipse] keys on the mesh, docs/geometry/Shapes.md"
 * The thickness is centred on the centre cell. The outline is a polygon of 360 vertices.
 * Staircasing: all three Yee points inside.
 *
 * @param radius_x, radius_z  semi-axes before the rotation [m]
 * @param thickness           along y [m]
 * @param theta               rotation in the x-z plane [rad]
 * @param antenna_shift       shift towards -y [m], whole cells
 */
bool fdtd::geometry::create_ellipse(int centerx, int centery_bulk, int centerz, double radius_x, double radius_z, double thickness, double theta, int index, double antenna_shift)
{
	if (!sc.draw || !(radius_x > 0 && radius_z > 0))
		return true;
	std::vector<std::pair<double, double>> outline;
	for (int n = 0; n < 360; ++n)
		outline.push_back({radius_x * std::cos(n * PI / 180), radius_z * std::sin(n * PI / 180)});
	const geom::Shape ellipse = geom::polygon(outline).rotated(theta, 0, 0);
	render_parts({{geom::extrude(ellipse, -0.5 * thickness, 0.5 * thickness), index}}, centerx, centery_bulk - int(antenna_shift / lat.d[1]), centerz);
	return true;
}

/*
*/
