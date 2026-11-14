#pragma once
/**
 * @file geometry/anapole.h
 * @brief Plasmonic anapole nanoantenna drawn with the geometry layers, [anapole], shape = anapole
 */
#include "geometry/csg.h"
#include <vector>

namespace geom
{
/// the plasmonic anapole: patch with a two-lobe aperture and a neck, [anapole] in the configuration
struct AnapoleSpec
{
	double patch_edge = 400e-9;          ///< [m]
	double thickness = 30e-9;            ///< [m] film thickness along y
	double corner_radius = 22e-9;        ///< [m] fillet of the patch corners
	double lobe_a = 0.433641, lobe_b = 0.445090, lobe_p = 1.4098, lobe_x0 = -0.120141; ///< in units of patch_edge
	double neck_half_length = 6.0216e-9; ///< [m] along the lobe axis
	double neck_half_height = 7.957e-9;  ///< [m] across the lobe axis
	double rotation = 0;                 ///< [deg] in the x-z plane
	int    lobe_material = -1;           ///< fill of the lobes; -1 = background
	int    neck_material = -1;           ///< fill of the neck; -1 = background
	int    dimer_material = -1;          ///< two cylinders in the neck; -1 = none
	double dimer_radius = 2.5e-9, dimer_gap = 2e-9; ///< [m]
	int    lobe_points = 512;            ///< polygon vertices per lobe edge
	bool operator==(const AnapoleSpec&) const = default;
};

/// one body and the material it is rendered with
struct Part
{
	Solid body;
	int material;
};

/// The anapole as parts in drawing order, later parts overwrite earlier ones: the patch minus the
/// aperture, then the lobe and neck fills and the dimer when they have a material.
/// (cx, cy, cz) is the centre of the patch [m]; the antenna lies in the x-z plane, thickness along y.
/// geometric cross section [m^2]: patch area minus aperture area, analytic for the patch and the
/// neck, shoelace area of the lobe polygons
double anapole_cross_section(const AnapoleSpec& a);

std::vector<Part> anapole_parts(const AnapoleSpec& a, int patch_material, double cx, double cy, double cz);
} // namespace geom
