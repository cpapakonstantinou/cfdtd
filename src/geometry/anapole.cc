/**
 * @file geometry/anapole.cc
 * @brief Plasmonic anapole nanoantenna, port of the COMSOL model plasmonic_anapole.mph
 *
 * Drawing, in the local frame, u along the lobe axis, v across it, origin at the patch centre:
 *
 *     patch    = rounded_rect(w, w, corner_radius)
 *     lobe     = polygon of the curve v = +-w B s^p sqrt(1 - s^2), s = (u / w - X0) / A,
 *                for neck_half_length <= u <= w (X0 + A), closed by the neck face u = neck_half_length
 *     aperture = lobe + mirror(lobe) + rect(2 neck_half_length, 2 neck_half_height)
 *     antenna  = extrude(patch - aperture, thickness), rotated by `rotation` in the x-z plane
 *
 * Optional fills: lobe and neck materials, and an ITO dimer, two cylinders of radius
 * dimer_radius, gap dimer_gap, centred in the neck along v.
 *
 * Axis mapping COMSOL -> FDTD: COMSOL x, lobe axis, -> u -> FDTD x at rotation 0;
 * COMSOL y, neck / dimer axis, -> v -> FDTD z; COMSOL z, thickness, -> FDTD y.
 * Default parameters: the COMSOL revision with s0 = 2, w = 400 nm, corner radius 22 nm,
 * t = 30 nm, fixed neck 12.04 x 15.91 nm.
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
	const double w = a.patch_edge, tx = a.neck_half_length;
	const double x_tip = w * (a.lobe_x0 + a.lobe_a);
	auto half = [&](double u) {
		const double s = (u / w - a.lobe_x0) / a.lobe_a;
		return (s <= 0 || s >= 1) ? 0.0 : w * a.lobe_b * std::pow(s, a.lobe_p) * std::sqrt(1 - s * s);
	};

	// one lobe (u > 0): upper edge from the neck face to the tip, lower edge back. The vertices are
	// spaced as sin() in s, dense near the tip where the edge sqrt(1 - s^2) turns vertical
	std::vector<std::pair<double, double>> pts;
	const int n = a.lobe_points;
	const double s_neck = (tx / w - a.lobe_x0) / a.lobe_a;
	auto u_at = [&](int m) {
		const double s = s_neck + (1 - s_neck) * std::sin(0.5 * M_PI * m / n);
		return m == n ? x_tip : w * (a.lobe_x0 + a.lobe_a * s);
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
	const double w = a.patch_edge, r = a.corner_radius;
	const double patch = w * w - (4 - M_PI) * r * r;
	const double neck = 4 * a.neck_half_length * a.neck_half_height;
	return patch - 2 * lobe_polygon(a, 0, 0).area() - neck; // lobes and neck meet only at the neck face
}

std::vector<Part> anapole_parts(const AnapoleSpec& a, int patch_material, double cx, double cy, double cz)
{
	const double w = a.patch_edge, tx = a.neck_half_length, nh = a.neck_half_height;
	const Shape lobe = lobe_polygon(a, cx, cz);
	const Shape lobes = lobe + lobe.mirrored_x(cx);
	const Shape neck = rect(cx, cz, 2 * tx, 2 * nh);
	const Shape patch = rounded_rect(cx, cz, w, w, a.corner_radius);

	const double theta = a.rotation * M_PI / 180.0;
	const double y0 = cy - 0.5 * a.thickness, y1 = cy + 0.5 * a.thickness;
	auto body = [&](const Shape& s) { return extrude(s.rotated(theta, cx, cz), y0, y1); };

	std::vector<Part> parts;
	parts.push_back({body(patch - (lobes + neck)), patch_material});
	if (a.lobe_material >= 0)
		parts.push_back({body(lobes - neck), a.lobe_material});
	if (a.neck_material >= 0)
		parts.push_back({body(neck), a.neck_material});
	if (a.dimer_material >= 0 && a.dimer_radius > 0)
	{
		const double off = 0.5 * a.dimer_gap + a.dimer_radius;
		parts.push_back({body(circle(cx, cz + off, a.dimer_radius) + circle(cx, cz - off, a.dimer_radius)), a.dimer_material});
	}
	return parts;
}

} // namespace geom
