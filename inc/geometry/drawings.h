#pragma once
/**
 * @file geometry/drawings.h
 * @brief The shapes drawn with the geometry layers, inc/geometry/csg.h: each a list of parts, a solid
 *        and its material, in drawing order; geometry/render.h writes them into the material maps
 *
 * Every drawing is in metres about the structure centre, the shape in the x-z plane and its
 * thickness along y unless it says otherwise. They need no solver state, test/test_geometry.cc.
 */
#include "geometry/csg.h"
#include <vector>

namespace geom
{
/// one body and the material it is rendered with
struct Part
{
	Solid body;
	int material;
};

/// shape = cylinder_shell: a core inside a shell, both along y; a length of 0 runs through the domain. The
/// core overwrites the shell where they overlap.
std::vector<Part> cylinder_shell_parts(double radius_core, double core_length, double radius_shell, double shell_length, int core_material, int shell_material);

/// shape = ellipse: an elliptic disc, semi-axes radius_x and radius_z turned by theta [rad], thickness centred on y = 0
std::vector<Part> ellipse_parts(double radius_x, double radius_z, double thickness, double theta, int material);

/// shape = anapole, after Hassan et al., Laser Photonics Rev. 19, 2400118, 2025: a patch with a dumbbell
/// aperture, two lobes joined by a gap, [anapole] in the configuration
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

/// The anapole as parts in drawing order, later parts overwrite earlier ones: the patch minus the
/// aperture, then the aperture fill and the gap fill when they have a material.
/// (cx, cy, cz) is the centre of the patch [m]; the antenna lies in the x-z plane, thickness along y.
/// geometric cross section [m^2]: patch area minus aperture area, analytic for the patch and the
/// gap, shoelace area of the lobe polygons
double anapole_cross_section(const AnapoleSpec& a);

std::vector<Part> anapole_parts(const AnapoleSpec& a, int patch_material, double cx, double cy, double cz);
} // namespace geom
