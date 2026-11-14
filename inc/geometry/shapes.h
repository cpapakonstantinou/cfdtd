/**
 * @file geometry/shapes.h
 * @brief The shape codes: the values of `[geometry] shape`, names as in the INI file
 */
#ifndef FDTD_GEOMETRY_SHAPES_H
#define FDTD_GEOMETRY_SHAPES_H

namespace fdtd::geometry::shape
{
inline constexpr int none = 0;
inline constexpr int sphere = 1;
inline constexpr int dipole = 5;
inline constexpr int bowtie = 6;
inline constexpr int cylinder = 12;
inline constexpr int two_spheres = 13;
inline constexpr int l_antenna = 17;
inline constexpr int hole = 19;
inline constexpr int lateral_dipole = 50;
inline constexpr int round_butterfly = 54;
inline constexpr int round_dipole = 55;
inline constexpr int one_side_round_dipole = 56;
inline constexpr int round_tripole = 57;
inline constexpr int round_hybrid = 58;
inline constexpr int x_shape = 59;
inline constexpr int sphere_dimer = 67;
inline constexpr int sphere_trimer = 68;
inline constexpr int sphere_quadrumer = 69;
inline constexpr int ellipse = 86;
inline constexpr int sphere_shell = 100;
inline constexpr int triangle = 103;
inline constexpr int rounded_triangle = 104;
inline constexpr int rounded_monopole = 105;
inline constexpr int u_shape = 106;
inline constexpr int cylinder_3d = 108;
inline constexpr int h_shape = 109;
inline constexpr int v_shape = 110;
inline constexpr int l_shape = 111;
inline constexpr int two_films_spacer = 112;
inline constexpr int rounded_dipole = 113;
inline constexpr int rounded_dipole_offset_planes = 114;
inline constexpr int cylinder_shell = 115;
inline constexpr int anapole = 2001;
} // namespace fdtd::geometry::shape

#endif
