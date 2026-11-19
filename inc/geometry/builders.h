/**
 * @file geometry/builders.h
 * @brief The shape builders: each writes its material slot into the material maps (matind_x/y/z)
 *        over the cells its shape covers, on this rank's part of the grid
 *
 * The arguments are the shape's centre in cells, its dimensions in m, the material slot and the layer
 * offsets of the run, in the order of the legacy builders; the set-arithmetic shapes of
 * inc/geometry.h, namespace geom, replace them one by one, roadmap §8. Returns false on failure.
 */
#ifndef FDTD_GEOMETRY_BUILDERS_H
#define FDTD_GEOMETRY_BUILDERS_H

namespace fdtd::geometry
{
/// basic.cc: primitives and the classic antennas
bool create_sphere(const int &, const int &, const int &, const double &, const int &, const double &);
bool create_sphere_2(const double &, const int &, const double &, const double &, const int &, const double &);
bool create_dipole(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const double &);
bool create_bowtie(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const int &, const double &);
bool create_rotated_monopole(const int &, const int &, const int &, const double &, const double &, const double &, const int &, const double &, const double &, const double &, const double &);
bool create_triangle(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const int &, const double &);
bool create_rounded_triangle(const int &, const int &, const int &, const double &, const double &, const double &, const int &, const double &, const double &, const double &);
double rounded_triangle_area(double w, double l, double r);
bool create_round_monopole(const int &, const int &, const int &, const double &, const double &, const double &, const int &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &);
bool create_monopole_rounded_lateral(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const int &, const double &);
bool create_round_dipole(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &);
bool create_layer(const int &, const int &, const int &);
bool create_anisotropic_layer(const int &, const int &, const int &, const int &, const int &);
bool create_cylinder(const int &, const int &, const int &, const double &, const int &, double length = 0);
/// the [block.*] boxes, in file order, blocks, set by apply_sim_config
bool create_cylinder_shell(int centerx, int centery, int centerz, double radius_core, double core_length, double radius_shell, double shell_length, int index, int shell_index);
bool create_ellipse(int centerx, int centery_bulk, int centerz, double radius_x, double radius_z, double thickness, double theta, int index, double antenna_shift);
bool create_cylinder_xz(const int &, const int &, const int &, const double &, const double &, const double &, const int &);

/// legacy.cc: older shapes kept on the keep-list
bool create_lateral_dipole(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const double &, const double &, const double &);
bool create_dipole_round(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const double &, const double &, const double &, const double &, const double &);
bool create_x_shaped(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const double &, const double &, const double &, const double &, const double &, const double &);
bool create_hole(const int &, const int &, const int &, const double &, const double &, const double &, const int &);

/// composite.cc: composite antennas
bool create_l_antenna(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const double &, const double &, const double &, const double &);
bool create_lh_butterfly(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const double &, const double &, const double &, const double &, const double &);
bool create_round_tripole(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const double &, const double &, const double &, const double &, const double &);
bool create_round_hybrid(const int &, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const double &, const double &, const double &, const double &, const double &);
bool create_background(const int &, const int &, const int &, const double &, const double &, const double &, const int &);

} // namespace fdtd::geometry

#endif
