/**
 * @file app/structures.cc
 * @brief The structure of the run from the configuration: which builders draw what, app/structures.h
 *
 * Application policy, like the plane placement, app/planes.cc: the shape named in [geometry] and its
 * dimensions are mapped onto calls of the geometry builders, inc/geometry/builders.h, after the
 * background layers and the [block.*] sections. Moved out of the geometry module, was
 * geometry/read_in.cc, on 2026-11-13.
 */
#include <stdio.h>
#include <cstring>
#include "materials/material_map.h"
#include "sources/signal.h"
#include "grid/subdomain.h"
#include "grid/state.h"
#include "app/run_config.h"
#include "monitors/state.h"
#include "materials/nonlinear/state.h"
#include "sources/state.h"
#include "geometry/state.h"
#include "io/text.h"
#include "geometry/builders.h"
#include "materials/materials.h"
#include "geometry/render.h"
#include "app/structures.h"
#include "geometry/scene.h"
#include "io/binary.h"
#include "grid/lattice.h"
using fdtd::constants::C0;
using fdtd::constants::PI;
using fdtd::constants::PI2;
using fdtd::grid::subdomain;
using fdtd::materials::matind_x;
using fdtd::materials::matind_y;
using fdtd::materials::matind_z;

namespace
{
const fdtd::app::RunConfig &run_cfg = fdtd::app::run_config(); ///< the run configuration
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::geometry; // the builders, the shape codes and the geometry state

geom::AnapoleSpec fdtd::app::anapole_input;
std::vector<fdtd::geometry::Block> fdtd::app::blocks;

void fdtd::app::create_structures()
{
	{
		fdtd::geometry::Scene sc;
		sc.shape = run_cfg.geometry.antenna_type;
		sc.draw = run_cfg.geometry.nano_antenna;
		sc.staircase = run_cfg.geometry.sta_ir_casing;
		sc.medium[0] = run_cfg.materials.first_medium, sc.medium[1] = run_cfg.materials.middle_medium, sc.medium[2] = run_cfg.materials.second_medium;
		sc.centre[0] = run_cfg.geometry.center_position_x, sc.centre[1] = run_cfg.geometry.center_position_y, sc.centre[2] = run_cfg.geometry.center_position_z;
		sc.gap_material = 0; // the round_hybrid gap is vacuum
		fdtd::geometry::set_scene(sc);
	}
	const int code = run_cfg.materials.antenna_material;
	fdtd::monitors::geometric_cross_section = -1; // declared by the shape drawn below

	fdtd::geometry::create_background(run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_z, run_cfg.geometry.middle_layer, run_cfg.geometry.thickness, run_cfg.geometry.antenna_shift, run_cfg.materials.antenna_material);
	create_blocks(blocks);
	MPI_Barrier(MPI_COMM_WORLD);

	if (run_cfg.geometry.is(shape::sphere_shell) || (add_free_electron_layer && run_cfg.geometry.is(shape::sphere)))
	{

		fdtd::geometry::create_sphere(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width + lat.d[0] * 13, run_cfg.geometry.shell_material, run_cfg.geometry.antenna_shift + lat.d[0] * 5);
		fdtd::geometry::create_sphere(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.materials.antenna_material, run_cfg.geometry.antenna_shift + lat.d[0] * 5);
	}

	if (run_cfg.geometry.is(shape::round_butterfly))
		fdtd::geometry::create_lh_butterfly(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, antenna_gap, run_cfg.geometry.length_2, run_cfg.materials.antenna_material, run_cfg.geometry.middle_layer, run_cfg.geometry.antenna_shift, run_cfg.geometry.length_1_l, run_cfg.geometry.length_2_b, run_cfg.geometry.rotation_deg * PI / 180); //round everywhere

	else if (run_cfg.geometry.is(shape::l_antenna))
		fdtd::geometry::create_l_antenna(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, antenna_gap, run_cfg.geometry.length_2, run_cfg.materials.antenna_material, run_cfg.geometry.middle_layer, run_cfg.geometry.antenna_shift, run_cfg.geometry.length_1_l, run_cfg.geometry.length_2_b);

	else if (run_cfg.geometry.is(shape::sphere))
	{
		fdtd::geometry::create_sphere(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, code, run_cfg.geometry.antenna_shift);
		fdtd::monitors::geometric_cross_section = PI * run_cfg.geometry.width * run_cfg.geometry.width;
	}
	else if (run_cfg.geometry.is(shape::cylinder_3d))
	{
		fdtd::geometry::create_cylinder_xz(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.length, run_cfg.geometry.antenna_shift, run_cfg.materials.antenna_material);
		fdtd::monitors::geometric_cross_section = PI * run_cfg.geometry.width * run_cfg.geometry.width;
	}
	else if (run_cfg.geometry.is(shape::dipole))
	{
		fdtd::geometry::create_dipole(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, antenna_gap, run_cfg.geometry.length_2, run_cfg.materials.antenna_material, run_cfg.geometry.antenna_shift);
		fdtd::monitors::geometric_cross_section = run_cfg.geometry.width * (run_cfg.geometry.length + run_cfg.geometry.length_2 - antenna_gap); // two arms: [g/2, l_t] and [-l_b, -g/2]
	}

	else if (run_cfg.geometry.is(shape::round_dipole))
	{
		fdtd::geometry::create_round_dipole(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, antenna_gap, run_cfg.geometry.length_2, run_cfg.materials.antenna_material, run_cfg.geometry.middle_layer, run_cfg.geometry.antenna_shift, run_cfg.geometry.length_1_l, run_cfg.geometry.length_2_b, run_cfg.geometry.rotation_deg * PI / 180, corner_radius, corner_radius, corner_radius); // rounded everywhere
		// two rounded blocks seen from above: w x l minus the four corner roundings
		auto arm = [](double w, double l, double r) { return (w > 0 && l > 0) ? w * l - (4 - PI) * r * r : double(0); };
		fdtd::monitors::geometric_cross_section = arm(run_cfg.geometry.width, run_cfg.geometry.length_2, corner_radius) + arm(run_cfg.geometry.width, run_cfg.geometry.length_2_b, corner_radius);
	}
	// bool fdtd::geometry::create_round_dipole(const int &centerx,const int &centery_bulk,const int &centerz,
	// const double &w, const double &t, const double &l1_r, const double &g,const double &l2_t,const int &index, const double &middle_layer, const double &antenna_shift,const double &l1_l, const double &l2_b, const double &alpha, const double &round_centre, const double &round_near, const double &round_far)
	else if (run_cfg.geometry.is(shape::round_tripole))
		fdtd::geometry::create_round_tripole(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, antenna_gap, run_cfg.geometry.length_2, run_cfg.materials.antenna_material, run_cfg.geometry.middle_layer, run_cfg.geometry.antenna_shift, run_cfg.geometry.length_1_l, run_cfg.geometry.length_2_b, run_cfg.geometry.rotation_deg * PI / 180); // rounded everywhere

	else if (run_cfg.geometry.is(shape::round_hybrid))
	{
		fdtd::geometry::create_round_hybrid(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, antenna_gap, run_cfg.geometry.length_2, code, run_cfg.geometry.middle_layer, run_cfg.geometry.antenna_shift, run_cfg.geometry.length_1_l, run_cfg.geometry.length_2_b, run_cfg.geometry.rotation_deg * PI / 180); // rounded everywhere
	}

	else if (run_cfg.geometry.is(shape::one_side_round_dipole))
	{
		fdtd::geometry::create_dipole_round(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, antenna_gap, run_cfg.geometry.length_2, code, run_cfg.geometry.middle_layer, run_cfg.geometry.antenna_shift, run_cfg.geometry.length_1_l, run_cfg.geometry.length_2_b, run_cfg.geometry.rotation_deg * PI / 180); // only extreme rounded
	}

	else if (run_cfg.geometry.is(shape::x_shape))
		fdtd::geometry::create_x_shaped(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, antenna_gap, run_cfg.geometry.length_2, run_cfg.materials.antenna_material, run_cfg.geometry.middle_layer, run_cfg.geometry.antenna_shift, run_cfg.geometry.length_1_l, run_cfg.geometry.length_2_b, run_cfg.geometry.rotation_deg * PI / 180, corner_radius);

	else if (run_cfg.geometry.is(shape::bowtie))
	{
		fdtd::geometry::create_bowtie(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, antenna_gap, code, run_cfg.geometry.antenna_shift);
		fdtd::monitors::geometric_cross_section = 0.5 * run_cfg.geometry.length * run_cfg.geometry.width - 0.5 * antenna_gap * antenna_gap * (run_cfg.geometry.width / run_cfg.geometry.length);
	}
	else if (run_cfg.geometry.is(shape::triangle))
	{
		fdtd::geometry::create_triangle(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, antenna_gap, code, run_cfg.geometry.antenna_shift);
		fdtd::monitors::geometric_cross_section = 0.5 * run_cfg.geometry.width * run_cfg.geometry.length;
	}
	else if (run_cfg.geometry.is(shape::rounded_triangle))
	{
		const double r = run_cfg.geometry.length_1_l; // corner radius; width x length is the rounded outline
		fdtd::geometry::create_rounded_triangle(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, code, run_cfg.geometry.antenna_shift, r, 0.0);
		fdtd::monitors::geometric_cross_section = fdtd::geometry::rounded_triangle_area(run_cfg.geometry.width, run_cfg.geometry.length, r);
	}
	else if (run_cfg.geometry.is(shape::rounded_monopole) || (run_cfg.geometry.is(shape::rounded_dipole) || run_cfg.geometry.is(shape::rounded_dipole_offset_planes)))
	{
		// rounded_monopole: one w x l block; rounded_dipole: two, separated by the gap along z
		const double r = run_cfg.geometry.length_1_l;
		const int arm = int((antenna_gap * 0.5 + run_cfg.geometry.length * 0.5) / lat.d[2]);
		for (int side : (run_cfg.geometry.is(shape::rounded_dipole) || run_cfg.geometry.is(shape::rounded_dipole_offset_planes)) ? std::vector<int>{1, -1} : std::vector<int>{0})
			fdtd::geometry::create_monopole_rounded_lateral(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z + side * arm, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, r, code, run_cfg.geometry.antenna_shift);
		fdtd::monitors::geometric_cross_section = ((run_cfg.geometry.is(shape::rounded_dipole) || run_cfg.geometry.is(shape::rounded_dipole_offset_planes)) ? 2 : 1) * (run_cfg.geometry.width * run_cfg.geometry.length - (4 - PI) * r * r);
	}

	else if (run_cfg.geometry.is(shape::u_shape))
	{
		double POINTRADIUS;
		POINTRADIUS = run_cfg.geometry.length_1_l;

		if (POINTRADIUS > run_cfg.geometry.width / 2.0)
			POINTRADIUS = run_cfg.geometry.width / 2.0;
		if (POINTRADIUS > run_cfg.geometry.thickness / 2.0)
			POINTRADIUS = run_cfg.geometry.thickness / 2.0;

		fdtd::geometry::create_round_monopole(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y + (run_cfg.geometry.thickness / 2.0) / lat.d[1], run_cfg.geometry.center_position_z - (run_cfg.geometry.length / 2.0 - run_cfg.geometry.width / 2.0) / lat.d[2], run_cfg.geometry.length_2 + run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.width, code, 0, run_cfg.geometry.antenna_shift, PI / 4, 0, POINTRADIUS, 0, 0, 0, 0);
		fdtd::geometry::create_round_monopole(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y + (run_cfg.geometry.thickness / 2.0) / lat.d[1], run_cfg.geometry.center_position_z + (run_cfg.geometry.length / 2.0 - run_cfg.geometry.width / 2.0) / lat.d[2], run_cfg.geometry.length_2 + run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.width, code, 0, run_cfg.geometry.antenna_shift, PI / 4, 0, POINTRADIUS, 0, 0, 0, 0);
		fdtd::geometry::create_round_monopole(run_cfg.geometry.center_position_x - (run_cfg.geometry.length_2 / lat.d[0]) / 2.0, run_cfg.geometry.center_position_y + (run_cfg.geometry.thickness / 2.0) / lat.d[1], run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length - run_cfg.geometry.width + POINTRADIUS * 2, code, 0, run_cfg.geometry.antenna_shift, PI / 4.0, 0, POINTRADIUS, 0, 0, 0, 0);
	}
	else if (run_cfg.geometry.is(shape::h_shape))
	{
		double POINTRADIUS;
		double DELTA;
		POINTRADIUS = run_cfg.geometry.length_1_l;
		DELTA = (run_cfg.geometry.length_2 / 2 - run_cfg.geometry.width / 2) / lat.d[0];
		// centery_bulk-ceil(middle_layer/DY)

		fdtd::geometry::create_layer(run_cfg.geometry.center_position_y - ceil((run_cfg.geometry.middle_layer + 200e-9) / lat.d[1]), run_cfg.geometry.center_position_y - ceil(run_cfg.geometry.middle_layer / lat.d[1]), run_cfg.materials.second_medium);

		fdtd::geometry::create_round_monopole(run_cfg.geometry.center_position_x + int(DELTA), run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, run_cfg.materials.antenna_material, 0, run_cfg.geometry.antenna_shift, PI / 4.0, 0, POINTRADIUS, 0, 0, 0, 0);
		fdtd::geometry::create_round_monopole(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.length_2 - run_cfg.geometry.width, run_cfg.geometry.thickness, antenna_gap, run_cfg.materials.antenna_material, 0, run_cfg.geometry.antenna_shift, PI / 4.0, 0, POINTRADIUS, 0, 0, 0, 0);
		fdtd::geometry::create_round_monopole(run_cfg.geometry.center_position_x - int(DELTA), run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, run_cfg.materials.antenna_material, 0, run_cfg.geometry.antenna_shift, PI / 4.0, 0, POINTRADIUS, 0, 0, 0, 0);

	}
	else if (run_cfg.geometry.is(shape::v_shape))
	{
		double POINTRADIUS;
		double DELTA;
		double ANGLE;
		ANGLE = PI / 4.0;

		POINTRADIUS = run_cfg.geometry.length_1_l;

		DELTA = ((run_cfg.geometry.width / 2.0 + run_cfg.geometry.length / 2.0) * cos(ANGLE) - run_cfg.geometry.width / (2.0 * sin(ANGLE))) / lat.d[0];
		// fdtd::io::log("geometry", "%e",Delta)
		// centery_bulk-ceil(middle_layer/DY)

		fdtd::geometry::create_rotated_monopole(run_cfg.geometry.center_position_x + int(DELTA) + 1, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, run_cfg.materials.antenna_material, 0, run_cfg.geometry.antenna_shift, ANGLE, POINTRADIUS);
		fdtd::geometry::create_rotated_monopole(run_cfg.geometry.center_position_x - int(DELTA) - 1, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, run_cfg.materials.antenna_material, 0, run_cfg.geometry.antenna_shift, ANGLE + PI / 2.0, POINTRADIUS);
		fdtd::geometry::create_rotated_monopole(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z - ((run_cfg.geometry.width / 2.0 + run_cfg.geometry.length / 2.0) / lat.d[2]) * cos(ANGLE) + sqrt(2) * run_cfg.geometry.width / lat.d[2], run_cfg.geometry.width, run_cfg.geometry.thickness, 2 * POINTRADIUS * cos(ANGLE), run_cfg.materials.antenna_material, 0, run_cfg.geometry.antenna_shift, 0, 0);
		fdtd::geometry::create_cylinder_xz(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y - ceil(run_cfg.geometry.antenna_shift / lat.d[1]), run_cfg.geometry.center_position_z - ((run_cfg.geometry.width / 2.0 + run_cfg.geometry.length / 2.0) / lat.d[2]) * cos(ANGLE) + sqrt(2) * run_cfg.geometry.width / lat.d[2] + 2 * (POINTRADIUS / lat.d[2]) * cos(ANGLE), POINTRADIUS, run_cfg.geometry.thickness, 0, run_cfg.materials.middle_medium);

	}
	else if (run_cfg.geometry.is(shape::l_shape))
	{
		double POINTRADIUS;
		POINTRADIUS = run_cfg.geometry.length_1_l;

		fdtd::geometry::create_rotated_monopole(run_cfg.geometry.center_position_x - int(run_cfg.geometry.length / (2.0 * lat.d[0])) + int(run_cfg.geometry.width / (2.0 * lat.d[0])), run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, run_cfg.materials.antenna_material, 0, run_cfg.geometry.antenna_shift, 0.0, POINTRADIUS);
		fdtd::geometry::create_rotated_monopole(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z - int(run_cfg.geometry.length / (2.0 * lat.d[0])) + int(run_cfg.geometry.width / (2.0 * lat.d[0])), run_cfg.geometry.length, run_cfg.geometry.thickness, run_cfg.geometry.width, run_cfg.materials.antenna_material, 0, run_cfg.geometry.antenna_shift, 0.0, POINTRADIUS);

	}

	else if (run_cfg.geometry.is(shape::cylinder))
	{
		fdtd::geometry::create_cylinder(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, code, run_cfg.geometry.length);
	}
	else if (run_cfg.geometry.is(shape::cylinder_shell))
	{
		// [cylinder_shell] radius_core, core_length, radius_shell, shell_length
		fdtd::geometry::create_cylinder_shell(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.length,
											  run_cfg.geometry.length_2, run_cfg.geometry.length_1_l, code, run_cfg.geometry.shell_material);
	}
	else if (run_cfg.geometry.is(shape::ellipse))
	{
		// [ellipse] radius_x, radius_z, thickness, rotation
		fdtd::geometry::create_ellipse(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.length,
									   run_cfg.geometry.thickness, run_cfg.geometry.rotation_deg * PI / 180, code, run_cfg.geometry.antenna_shift);
		fdtd::monitors::geometric_cross_section = PI * run_cfg.geometry.width * run_cfg.geometry.length;
	}
	else if (run_cfg.geometry.is(shape::hole))
	{
		fdtd::geometry::create_hole(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.antenna_shift, code);
	}

	else if (run_cfg.geometry.is(shape::two_spheres))
	{
		fdtd::geometry::create_sphere(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y - int(run_cfg.geometry.antenna_shift / lat.d[1] / 2), run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.materials.antenna_material, 0);
		fdtd::geometry::create_sphere(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y + int(run_cfg.geometry.antenna_shift / lat.d[1] / 2), run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.materials.antenna_material, 0);
	}

	else if (run_cfg.geometry.is(shape::sphere_dimer))
	{
		double d_z = run_cfg.geometry.length; //interdistance centre-to-centre
		double r_1 = run_cfg.geometry.width;
		double embed = run_cfg.geometry.antenna_shift;
		fdtd::geometry::create_sphere_2(double(run_cfg.geometry.center_position_x * lat.d[0]), run_cfg.geometry.center_position_y, double(run_cfg.geometry.center_position_z * lat.d[2] + 0.5 * d_z), r_1, run_cfg.materials.antenna_material, embed);
		fdtd::geometry::create_sphere_2(double(run_cfg.geometry.center_position_x * lat.d[0]), run_cfg.geometry.center_position_y, double(run_cfg.geometry.center_position_z * lat.d[2] - 0.5 * d_z), r_1, run_cfg.materials.antenna_material, embed);
	}

	else if (run_cfg.geometry.is(shape::sphere_trimer))
	{
		double d_z = run_cfg.geometry.length; //interdistance centre-to-centre
		double r_1 = run_cfg.geometry.width;
		double embed = run_cfg.geometry.antenna_shift;
		fdtd::geometry::create_sphere_2(double(run_cfg.geometry.center_position_x * lat.d[0]), run_cfg.geometry.center_position_y, double(run_cfg.geometry.center_position_z * lat.d[2]), r_1, run_cfg.materials.antenna_material, embed);
		fdtd::geometry::create_sphere_2(double(run_cfg.geometry.center_position_x * lat.d[0]), run_cfg.geometry.center_position_y, double(run_cfg.geometry.center_position_z * lat.d[2] + d_z), r_1, run_cfg.materials.antenna_material, embed);
		fdtd::geometry::create_sphere_2(double(run_cfg.geometry.center_position_x * lat.d[0]), run_cfg.geometry.center_position_y, double(run_cfg.geometry.center_position_z * lat.d[2] - d_z), r_1, run_cfg.materials.antenna_material, embed);
	}

	else if (run_cfg.geometry.is(shape::sphere_quadrumer))
	{
		double d_z = run_cfg.geometry.length; //interdistance centre-to-centre
		double r_1 = run_cfg.geometry.width;
		double embed = run_cfg.geometry.antenna_shift;
		fdtd::geometry::create_sphere_2(double(run_cfg.geometry.center_position_x * lat.d[0]), run_cfg.geometry.center_position_y, double(run_cfg.geometry.center_position_z * lat.d[2] + 0.5 * d_z), r_1, run_cfg.materials.antenna_material, embed);
		fdtd::geometry::create_sphere_2(double(run_cfg.geometry.center_position_x * lat.d[0]), run_cfg.geometry.center_position_y, double(run_cfg.geometry.center_position_z * lat.d[2] - 0.5 * d_z), r_1, run_cfg.materials.antenna_material, embed);
		fdtd::geometry::create_sphere_2(double(run_cfg.geometry.center_position_x * lat.d[0]), run_cfg.geometry.center_position_y, double(run_cfg.geometry.center_position_z * lat.d[2] + 1.5 * d_z), r_1, run_cfg.materials.antenna_material, embed);
		fdtd::geometry::create_sphere_2(double(run_cfg.geometry.center_position_x * lat.d[0]), run_cfg.geometry.center_position_y, double(run_cfg.geometry.center_position_z * lat.d[2] - 1.5 * d_z), r_1, run_cfg.materials.antenna_material, embed);
	}

	else if (run_cfg.geometry.is(shape::lateral_dipole))
		fdtd::geometry::create_lateral_dipole(run_cfg.geometry.center_position_x, run_cfg.geometry.center_position_y, run_cfg.geometry.center_position_z, run_cfg.geometry.width, run_cfg.geometry.thickness, run_cfg.geometry.length, antenna_gap, run_cfg.geometry.length_2, run_cfg.materials.antenna_material, run_cfg.geometry.middle_layer, run_cfg.geometry.antenna_shift, run_cfg.geometry.length_1_l);

	else if (run_cfg.geometry.is(shape::two_films_spacer))
	{
		// Puts a spacer in the Middle Medium.
		fdtd::geometry::create_layer(int(run_cfg.geometry.antenna_shift / lat.d[1]) + (run_cfg.geometry.center_position_y - int(run_cfg.geometry.middle_layer / lat.d[1])), int((run_cfg.geometry.thickness + run_cfg.geometry.antenna_shift) / lat.d[1]) + (run_cfg.geometry.center_position_y - int(run_cfg.geometry.middle_layer / lat.d[1])), run_cfg.materials.antenna_material);
	}

	else if (run_cfg.geometry.is(shape::anapole))
	{
		fdtd::geometry::create_anapole(anapole_input, run_cfg.materials.antenna_material, run_cfg.geometry.center_position_x,
							   run_cfg.geometry.center_position_y - int(run_cfg.geometry.antenna_shift / lat.d[1]), run_cfg.geometry.center_position_z);
		fdtd::monitors::geometric_cross_section = geom::anapole_cross_section(anapole_input);
	}
}

