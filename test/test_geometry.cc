// Unit tests of the draw-then-render geometry, inc/geometry.h.
// The raster of every drawing is checked cell by cell against direct point sampling.
#include "geometry/csg.h"
#include "geometry/anapole.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <array>
#include <random>

static int failures = 0;
#define CHECK(c) do { if (!(c)) { std::printf("FAIL %s:%d %s\n", __FILE__, __LINE__, #c); ++failures; } } while (0)

using Pred = std::function<bool(double, double)>; // is (x, z) inside?

static bool in_polygon(const std::vector<std::pair<double, double>>& p, double x, double z)
{
	// the rasteriser's rule: crossings of the vertical line at x, half-open in x
	std::vector<double> zs;
	for (size_t e = 0; e < p.size(); ++e)
	{
		const auto& a = p[e];
		const auto& b = p[(e + 1) % p.size()];
		if ((a.first <= x && x < b.first) || (b.first <= x && x < a.first))
			zs.push_back(a.second + (x - a.first) * (b.second - a.second) / (b.first - a.first));
	}
	std::sort(zs.begin(), zs.end());
	for (size_t c = 0; c + 1 < zs.size(); c += 2)
		if (zs[c] <= z && z <= zs[c + 1])
			return true;
	return false;
}

// raster vs. sampling the predicate at every lattice point
static long compare(const geom::Shape& s, const Pred& inside, const geom::Lattice& L)
{
	const geom::Layer layer = s.render(L);
	long bad = 0;
	for (int i = L.i0; i < L.i1; ++i)
	{
		const auto& row = layer.rows[i - L.i0];
		for (int k = L.k0; k < L.k1; ++k)
		{
			bool in_raster = false;
			for (const auto& r : row)
				in_raster = in_raster || (k >= r.k0 && k < r.k1);
			bad += in_raster != inside(L.x(i), L.z(k));
		}
	}
	return bad;
}


// the original point-wise anapole predicate, src/anapoleGeometry.cc draft: material or -1
static int draft_anapole(double u, double v, const geom::AnapoleSpec& a, int index)
{
	const double au = std::fabs(u);
	const double av = std::fabs(v);
	const double w = a.patch_edge;
	const double hw = 0.5 * w;
	const double tx = 0.5 * a.gap_length;
	const double nh = 0.5 * a.gap_width;
	const double x_tip = w * (a.lobe_offset + a.lobe_length);
	if (au <= tx && av <= nh)
		return a.gap_material;
	if (au >= tx && au <= x_tip)
	{
		const double s = (au / w - a.lobe_offset) / a.lobe_length;
		if (s > 0 && s <= 1 && av <= w * a.lobe_height * std::pow(s, a.lobe_exponent) * std::sqrt(1 - s * s))
			return a.aperture_material;
	}
	if (au <= hw && av <= hw)
	{
		const double c = hw - a.corner_radius;
		if (a.corner_radius > 0 && au > c && av > c && (au - c) * (au - c) + (av - c) * (av - c) > a.corner_radius * a.corner_radius)
			return -1;
		return index;
	}
	return -1;
}

int main()
{
	using namespace geom;

	// run algebra
	CHECK((unite({{0, 3}, {5, 7}}, {{2, 6}}) == Row{{0, 7}}));
	CHECK((unite({{0, 2}}, {{2, 4}}) == Row{{0, 4}}));
	CHECK((subtract({{0, 10}}, {{2, 3}, {5, 7}}) == Row{{0, 2}, {3, 5}, {7, 10}}));
	CHECK((subtract({{0, 4}}, {{0, 4}}).empty()));
	CHECK((intersect({{0, 5}, {8, 12}}, {{3, 9}}) == Row{{3, 5}, {8, 9}}));
	{
		std::mt19937 rng(7);
		std::uniform_int_distribution<int> u(0, 60);
		for (int t = 0; t < 2000; ++t)
		{
			auto make = [&] {
				std::vector<bool> bits(64, false);
				for (int r = 0; r < 4; ++r)
				{
					int a = u(rng);
					int b = a + u(rng) % 8;
					for (int k = a; k < b && k < 64; ++k)
						bits[k] = true;
				}
				Row row;
				for (int k = 0; k < 64; ++k)
					if (bits[k] && (row.empty() || row.back().k1 != k))
						row.push_back({k, k + 1});
					else if (bits[k])
						row.back().k1 = k + 1;
				return std::make_pair(bits, row);
			};
			auto [ba, ra] = make();
			auto [bb, rb] = make();
			auto bits_of = [](const Row& row) {
				std::vector<bool> b(64, false);
				for (const auto& r : row)
					for (int k = r.k0; k < r.k1; ++k)
						b[k] = true;
				return b;
			};
			auto ok = [&](const Row& got, auto op) {
				auto g = bits_of(got);
				for (int k = 0; k < 64; ++k)
					if (g[k] != op(ba[k], bb[k]))
						return false;
				for (size_t m = 1; m < got.size(); ++m) // sorted, disjoint, not touching
					if (got[m].k0 <= got[m - 1].k1)
						return false;
				return true;
			};
			CHECK(ok(unite(ra, rb), [](bool a, bool b) { return a || b; }));
			CHECK(ok(subtract(ra, rb), [](bool a, bool b) { return a && !b; }));
			CHECK(ok(intersect(ra, rb), [](bool a, bool b) { return a && b; }));
		}
	}

	// primitives and set operations vs. point sampling
	Lattice L;
	L.dx = 1e-9, L.dz = 1e-9, L.i0 = -5, L.i1 = 75, L.k0 = -5, L.k1 = 75;
	for (double ox : {0.0, 0.5})
		for (double oz : {0.0, 0.5})
		{
			L.off_x = ox, L.off_z = oz;
			const double cx = 33.3e-9;
			const double cz = 31.7e-9;
			auto in_rect = [](double x, double z, double cx, double cz, double wx, double wz) {
				return std::fabs(x - cx) <= 0.5 * wx && std::fabs(z - cz) <= 0.5 * wz;
			};
			auto in_circle = [](double x, double z, double cx, double cz, double r) {
				return (x - cx) * (x - cx) + (z - cz) * (z - cz) <= r * r;
			};
			Shape arm1 = rect(cx, cz, 50.3e-9, 10.1e-9);
			Shape arm2 = rect(cx, cz, 10.1e-9, 50.3e-9);
			Shape hole = circle(cx, cz, 4.3e-9);
			Pred p_arm1 = [&](double x, double z) { return in_rect(x, z, cx, cz, 50.3e-9, 10.1e-9); };
			Pred p_arm2 = [&](double x, double z) { return in_rect(x, z, cx, cz, 10.1e-9, 50.3e-9); };
			Pred p_hole = [&](double x, double z) { return in_circle(x, z, cx, cz, 4.3e-9); };
			CHECK(compare(arm1, p_arm1, L) == 0);
			CHECK(compare(hole, p_hole, L) == 0);
			CHECK(compare((arm1 + arm2) - hole, [&](double x, double z) { return (p_arm1(x, z) || p_arm2(x, z)) && !p_hole(x, z); }, L) == 0);
			CHECK(compare(arm1 & arm2, [&](double x, double z) { return p_arm1(x, z) && p_arm2(x, z); }, L) == 0);

			std::vector<std::pair<double, double>> tri = {{5.2e-9, 6.1e-9}, {61.7e-9, 12.4e-9}, {22.9e-9, 58.3e-9}};
			CHECK(compare(polygon(tri), [&](double x, double z) { return in_polygon(tri, x, z); }, L) == 0);
			CHECK(compare(polygon(tri) - hole, [&](double x, double z) { return in_polygon(tri, x, z) && !p_hole(x, z); }, L) == 0);
		}

	// areas
	CHECK(std::fabs(rect(0, 0, 2, 3).area() - 6) < 1e-12);
	CHECK(std::fabs(polygon({{0, 0}, {4, 0}, {0, 3}}).area() - 6) < 1e-12);
	CHECK((rect(0, 0, 1, 1) + circle(0, 0, 1)).area() < 0); // composite: mesh area only
	{
		// a fine raster converges to the analytic area
		Lattice F;
		F.dx = F.dz = 0.01e-9, F.off_x = F.off_z = 0.5, F.i0 = -3000, F.i1 = 3000, F.k0 = -3000, F.k1 = 3000;
		const double r = 25.1e-9;
		const double mesh = circle(0, 0, r).render(F).cells() * F.dx * F.dz;
		CHECK(std::fabs(mesh / (M_PI * r * r) - 1) < 1e-3);
		const Shape rr = rounded_rect(0, 0, 40e-9, 30e-9, 5e-9);
		const double exact = 40e-9 * 30e-9 - (4 - M_PI) * 5e-9 * 5e-9;
		CHECK(std::fabs(rr.render(F).cells() * F.dx * F.dz / exact - 1) < 1e-3);
	}

	// transforms
	{
		const Shape r = rect(10e-9, 0, 4e-9, 2e-9);
		const Shape q = r.rotated(M_PI / 2, 0, 0); // -> centred at (0, 10 nm), 2 x 4 nm
		Lattice M;
		M.dx = M.dz = 0.5e-9, M.i0 = -40, M.i1 = 40, M.k0 = -40, M.k1 = 40, M.off_x = M.off_z = 0.5;
		CHECK(compare(q, [](double x, double z) { return std::fabs(x) <= 1e-9 && std::fabs(z - 10e-9) <= 2e-9; }, M) == 0);
		CHECK(compare(r.mirrored_x(0), [](double x, double z) { return std::fabs(x + 10e-9) <= 2e-9 && std::fabs(z) <= 1e-9; }, M) == 0);
	}

	// solids: slices
	{
		const Solid ball = sphere(0, 0, 0, 10.3e-9);
		CHECK(ball.slice(11e-9).empty());
		CHECK(std::fabs(ball.slice(6e-9).area() - M_PI * (10.3e-9 * 10.3e-9 - 36e-18)) < 1e-24);
		const Solid film = extrude(rect(0, 0, 100.6e-9, 100.6e-9), -15e-9, 15e-9);
		CHECK(!film.slice(15e-9).empty() && film.slice(15.1e-9).empty());
		const Solid holey = film - ball;
		Lattice M;
		M.dx = M.dz = 1e-9, M.i0 = -60, M.i1 = 60, M.k0 = -60, M.k1 = 60;
		CHECK(compare(holey.slice(0), [](double x, double z) {
			return std::fabs(x) <= 50.3e-9 && std::fabs(z) <= 50.3e-9 && x * x + z * z > 10.3e-9 * 10.3e-9;
		}, M) == 0);
		// degenerate inputs are empty, not inverted or point-like
		CHECK(extrude(rect(0, 0, 1e-9, 1e-9), 2e-9, -2e-9).empty());
		CHECK(sphere(0, 0, 0, 0).empty() && cylinder_x(0, 0, 0, -1, 1).empty() && cylinder_z(0, 0, -1e-9).empty());
	}

	// finite cylinders, moved solids
	{
		Lattice M;
		M.dx = M.dz = 1e-9, M.i0 = -60, M.i1 = 60, M.k0 = -60, M.k1 = 60;
		const double r = 5.3e-9;
		const double y = 1.3e-9;
		const double cy = 0.2e-9;
		const double h = std::sqrt(r * r - (y - cy) * (y - cy));
		CHECK(compare(cylinder_x(cy, 3.1e-9, r, -20.5e-9, 17.3e-9).slice(y), [&](double x, double z) {
			return x >= -20.5e-9 && x <= 17.3e-9 && std::fabs(z - 3.1e-9) <= h;
		}, M) == 0);
		CHECK(compare(cylinder_z(-4.1e-9, cy, r, -11.5e-9, 23.3e-9).slice(y), [&](double x, double z) {
			return z >= -11.5e-9 && z <= 23.3e-9 && std::fabs(x + 4.1e-9) <= h;
		}, M) == 0);
		CHECK(cylinder_x(cy, 0, r, -1e-9, 1e-9).slice(cy + 5.4e-9).empty());
		// a ball moved by (5.3, 1.1, -7.7) nm, then turned by 0.6 rad about (2, 3) nm
		const double R = 9.7e-9;
		const double th = 0.6;
		const double ox = 2e-9;
		const double oz = 3e-9;
		const double ys = 2.0e-9;
		const Solid moved = sphere(0, 0, 0, R).translated(5.3e-9, 1.1e-9, -7.7e-9).rotated(th, ox, oz);
		CHECK(compare(moved.slice(ys), [&](double x, double z) {
			const double c = std::cos(-th), s = std::sin(-th), u = ox + c * (x - ox) - s * (z - oz), v = oz + s * (x - ox) + c * (z - oz);
			return (u - 5.3e-9) * (u - 5.3e-9) + (v + 7.7e-9) * (v + 7.7e-9) <= R * R - (ys - 1.1e-9) * (ys - 1.1e-9);
		}, M) <= 2); // rotation round-off on at most a couple of boundary samples
	}

	// rounded polygon: the set of points within r of the inner, offset, polygon
	{
		Lattice M;
		M.dx = M.dz = 0.5e-9, M.i0 = -80, M.i1 = 80, M.k0 = -80, M.k1 = 80;
		const double w = 40.3e-9;
		const double l = 50.7e-9;
		const double r = 4.1e-9;
		const std::vector<std::pair<double, double>> tri{{-0.5 * w, -0.5 * l}, {0.5 * w, -0.5 * l}, {0, 0.5 * l}};
		// inner triangle: the three edges moved inward by r, intersected pairwise
		std::vector<std::pair<double, double>> inner;
		auto line = [&](int e) { // n . p = d with n the inward unit normal, shifted by r
			const auto [ax, az] = tri[e];
			const auto [bx, bz] = tri[(e + 1) % 3];
			double nx = -(bz - az);
			double nz = bx - ax; // left normal: inward for a counter-clockwise polygon
			const double nl = std::hypot(nx, nz);
			nx /= nl, nz /= nl;
			return std::array<double, 3>{nx, nz, nx * ax + nz * az + r};
		};
		for (int v = 0; v < 3; ++v)
		{
			const auto p = line((v + 2) % 3);
			const auto q = line(v);
			const double det = p[0] * q[1] - p[1] * q[0];
			inner.push_back({(p[2] * q[1] - p[1] * q[2]) / det, (p[0] * q[2] - p[2] * q[0]) / det});
		}
		auto dist_seg = [](double x, double z, std::pair<double, double> a, std::pair<double, double> b) {
			const double ux = b.first - a.first;
			const double uz = b.second - a.second;
			const double t = std::clamp(((x - a.first) * ux + (z - a.second) * uz) / (ux * ux + uz * uz), 0.0, 1.0);
			return std::hypot(x - a.first - t * ux, z - a.second - t * uz);
		};
		auto rounded = [&](double x, double z) {
			if (in_polygon(inner, x, z))
				return true;
			double d = 1e9;
			for (int e = 0; e < 3; ++e)
				d = std::min(d, dist_seg(x, z, inner[e], inner[(e + 1) % 3]));
			return d <= r;
		};
		CHECK(compare(rounded_polygon(tri, r), rounded, M) == 0);
		CHECK(compare(rounded_polygon(tri, 0), [&](double x, double z) { return in_polygon(tri, x, z); }, M) == 0);
		// clockwise order gives the same shape
		CHECK(compare(rounded_polygon({tri[2], tri[1], tri[0]}, r), rounded, M) == 0);
		// growing the inner triangle by r gives the same set
		CHECK(compare(inflated_polygon(inner, r), rounded, M) == 0);
		CHECK(compare(inflated_polygon({inner[2], inner[1], inner[0]}, r), rounded, M) == 0);
	}


	// anapole: drawing vs. the original point-wise predicate
	{
		geom::AnapoleSpec a;
		CHECK(a.patch_edge == 200e-9 && a.gap_length == 12.0432e-9 && a.gap_width == 15.914e-9);
		// the 400 nm variant, corner radius scaled with the edge
		a.patch_edge = 400e-9;
		a.corner_radius = 22e-9;
		const double cx = 300.3e-9;
		const double cz = 299.7e-9;
		const auto parts = anapole_parts(a, 7, cx, 0, cz);
		CHECK(parts.size() == 1);
		Lattice M;
		M.dx = M.dz = 1e-9, M.off_x = 0.5, M.i0 = 80, M.i1 = 520, M.k0 = 80, M.k1 = 520;
		const Layer layer = parts[0].body.slice(0).render(M);
		long bad = 0;
		long ag = 0;
		for (int i = M.i0; i < M.i1; ++i)
			for (int k = M.k0; k < M.k1; ++k)
			{
				bool in = false;
				for (const auto& r : layer.rows[i - M.i0])
					in = in || (k >= r.k0 && k < r.k1);
				const bool want = draft_anapole(M.x(i) - cx, M.z(k) - cz, a, 7) == 7;
				bad += in != want;
				ag += want;
			}
		std::printf("anapole: %ld patch cells, %ld differ from the point-wise predicate\n", ag, bad);
		// the declared cross section matches the cell count of the 1 nm raster
		CHECK(std::fabs(anapole_cross_section(a) / (ag * M.dx * M.dz) - 1) < 2e-3);
		CHECK(ag > 100000 && bad <= 4); // lobe edge sampled with 512 vertices: chord error << 1 nm
		CHECK(parts[0].body.slice(16e-9).empty() && !parts[0].body.slice(15e-9).empty());
		a.aperture_material = 3;
		CHECK(anapole_parts(a, 7, cx, 0, cz).size() == 2);
		a.gap_material = 4;
		CHECK(anapole_parts(a, 7, cx, 0, cz).size() == 3);
	}

	if (failures == 0)
		std::printf("all geometry tests passed\n");
	return failures ? 1 : 0;
}
