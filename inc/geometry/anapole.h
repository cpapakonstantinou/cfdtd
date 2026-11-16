#pragma once
/**
 * @file geometry/anapole.h
 * @brief Plasmonic anapole nanoantenna drawn with the geometry layers, [anapole], shape = anapole
 */
#include "geometry/csg.h"
#include <vector>

namespace geom
{
/// the plasmonic anapole: patch with a dumbbell aperture, two lobes joined by a gap, [anapole] in the configuration
struct AnapoleSpec
{
	double patch_edge = 200e-9; ///< [m] edge w of the square patch
	double thickness = 30e-9; ///< [m] film thickness along y
	double corner_radius = 11e-9; ///< [m] fillet of the patch corners
	double lobe_length = 0.433641; ///< A, in units of patch_edge
	double lobe_height = 0.445090; ///< B, in units of patch_edge
	double lobe_exponent = 1.4098; ///< p
	double lobe_offset = -0.120141; ///< X0, where the lobe curve starts, in units of patch_edge
	double gap_length = 12.0432e-9; ///< [m] gap along the lobe axis
	double gap_width = 15.914e-9; ///< [m] gap across the lobe axis
	double rotation = 0; ///< [deg] in the x-z plane
	int aperture_material = -1; ///< fill of the whole aperture, lobes and gap; -1 = background
	int gap_material = -1; ///< fill of the gap alone, drawn over the aperture fill; -1 = as the aperture
	int lobe_points = 512; ///< polygon vertices per lobe edge
	bool operator==(const AnapoleSpec&) const = default;
};

/// one body and the material it is rendered with
struct Part
{
	Solid body;
	int material;
};

/// The anapole as parts in drawing order, later parts overwrite earlier ones: the patch minus the
/// aperture, then the aperture fill and the gap fill when they have a material.
/// (cx, cy, cz) is the centre of the patch [m]; the antenna lies in the x-z plane, thickness along y.
/// geometric cross section [m^2]: patch area minus aperture area, analytic for the patch and the
/// gap, shoelace area of the lobe polygons
double anapole_cross_section(const AnapoleSpec& a);

std::vector<Part> anapole_parts(const AnapoleSpec& a, int patch_material, double cx, double cy, double cz);
} // namespace geom
