/**
 * @file geometry/anapole.cc
 * @brief Plasmonic anapole nanoantenna after Hassan et al., Laser Photonics Rev. 19, 2400118, 2025
 *
 * Drawing, in the local frame, u along the lobe axis, v across it, origin at the patch centre:
 *
 *     patch    = rounded_rect(w, w, corner_radius)
 *     gap      = rect(gap_length, gap_width)
 *     lobe     = polygon of the curve v = +-w B s^p sqrt(1 - s^2), s = (|u| / w - X0) / A,
 *                from the gap face |u| = gap_length / 2 to the tip |u| = w (X0 + A)
 *     aperture = the two lobes + gap
 *     antenna  = extrude(patch - aperture, thickness), rotated by `rotation` in the x-z plane
 *
 * with w = patch_edge, A = lobe_length, B = lobe_height, p = lobe_exponent, X0 = lobe_offset.
 * Optional fills: aperture_material fills the whole aperture, gap_material the gap alone.
 *
 * At rotation 0 the lobe axis u is x and v is z; the thickness is along y.
 * Default parameters: the fit to the published design, w = 200 nm, t = 30 nm, gap 12.04 x 15.91 nm.
 */
#include "geometry/anapole.h"

#include <cmath>

namespace geom
{

namespace
{
/// one lobe (u > 0) about (cx, cz)
Shape lobe_polygon(const AnapoleSpec& a, double cx, double cz)
{
	const double w = a.patch_edge;
	const double tx = 0.5 * a.gap_length;
	const double x_tip = w * (a.lobe_offset + a.lobe_length);
	auto half = [&](double u) {
		const double s = (u / w - a.lobe_offset) / a.lobe_length;
		return (s <= 0 || s >= 1) ? 0.0 : w * a.lobe_height * std::pow(s, a.lobe_exponent) * std::sqrt(1 - s * s);
	};

	// one lobe (u > 0): upper edge from the gap face to the tip, lower edge back. The vertices are
	// spaced as sin() in s, dense near the tip where the edge sqrt(1 - s^2) turns vertical
	std::vector<std::pair<double, double>> pts;
	const int n = a.lobe_points;
	const double s_gap = (tx / w - a.lobe_offset) / a.lobe_length;
	auto u_at = [&](int m) {
		const double s = s_gap + (1 - s_gap) * std::sin(0.5 * M_PI * m / n);
		return m == n ? x_tip : w * (a.lobe_offset + a.lobe_length * s);
	};
	for (int m = 0; m <= n; ++m)
		pts.push_back({cx + u_at(m), cz + half(u_at(m))});
	for (int m = n; m >= 0; --m)
		pts.push_back({cx + u_at(m), cz - half(u_at(m))});
	return polygon(pts);
}
} // namespace

double anapole_cross_section(const AnapoleSpec& a)
{
	const double w = a.patch_edge;
	const double r = a.corner_radius;
	const double patch = w * w - (4 - M_PI) * r * r;
	const double gap = a.gap_length * a.gap_width;
	return patch - 2 * lobe_polygon(a, 0, 0).area() - gap; // lobes and gap meet only at the gap face
}

std::vector<Part> anapole_parts(const AnapoleSpec& a, int patch_material, double cx, double cy, double cz)
{
	const double w = a.patch_edge;
	const Shape lobe = lobe_polygon(a, cx, cz);
	const Shape lobes = lobe + lobe.mirrored_x(cx);
	const Shape gap = rect(cx, cz, a.gap_length, a.gap_width);
	const Shape patch = rounded_rect(cx, cz, w, w, a.corner_radius);

	const double theta = a.rotation * M_PI / 180.0;
	const double y0 = cy - 0.5 * a.thickness;
	const double y1 = cy + 0.5 * a.thickness;
	auto body = [&](const Shape& s) { return extrude(s.rotated(theta, cx, cz), y0, y1); };

	std::vector<Part> parts;
	parts.push_back({body(patch - (lobes + gap)), patch_material});
	if (a.aperture_material >= 0)
		parts.push_back({body(lobes + gap), a.aperture_material});
	if (a.gap_material >= 0)
		parts.push_back({body(gap), a.gap_material});
	return parts;
}

} // namespace geom
