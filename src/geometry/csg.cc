/**
 * @file geometry/csg.cc
 * @brief Draw-then-render geometry, see geometry/csg.h
 */
#include "geometry/csg.h"

#include <algorithm>
#include <cmath>

namespace geom
{

// runs

Row unite(const Row& a, const Row& b)
{
	Row out;
	out.reserve(a.size() + b.size());
	size_t i = 0, j = 0;
	while (i < a.size() || j < b.size())
	{
		const Run r = (j == b.size() || (i < a.size() && a[i].k0 <= b[j].k0)) ? a[i++] : b[j++];
		if (!out.empty() && r.k0 <= out.back().k1)
			out.back().k1 = std::max(out.back().k1, r.k1);
		else
			out.push_back(r);
	}
	return out;
}

Row subtract(const Row& a, const Row& b)
{
	Row out;
	size_t j = 0;
	for (Run r : a)
	{
		while (j < b.size() && b[j].k1 <= r.k0)
			++j;
		for (size_t m = j; m < b.size() && b[m].k0 < r.k1; ++m)
		{
			if (b[m].k0 > r.k0)
				out.push_back({r.k0, b[m].k0});
			r.k0 = std::max(r.k0, b[m].k1);
			if (r.k0 >= r.k1)
				break;
		}
		if (r.k0 < r.k1)
			out.push_back(r);
	}
	return out;
}

Row intersect(const Row& a, const Row& b)
{
	Row out;
	size_t i = 0, j = 0;
	while (i < a.size() && j < b.size())
	{
		const int lo = std::max(a[i].k0, b[j].k0), hi = std::min(a[i].k1, b[j].k1);
		if (lo < hi)
			out.push_back({lo, hi});
		(a[i].k1 < b[j].k1) ? ++i : ++j;
	}
	return out;
}

long Layer::cells() const
{
	long n = 0;
	for (const auto& row : rows)
		for (const auto& r : row)
			n += r.k1 - r.k0;
	return n;
}

// shapes

struct Shape::Node
{
	enum Kind { Rect, Circle, Polygon, Union, Difference, Intersection } kind = Rect;
	double x0 = 0, x1 = 0, z0 = 0, z1 = 0;            // rect: bounds
	double cx = 0, cz = 0, r = 0;                     // circle: centre, radius
	std::vector<std::pair<double, double>> pts;        // polygon
	Shape a, b;                                        // composites
};

namespace
{
const double EPS = 1e-9; // fraction of a cell: samples on a boundary are inside

/// the samples k with z_k in the closed interval [za, zb], clipped to the lattice
void add_interval(Row& row, double za, double zb, const Lattice& L)
{
	// clip in floating point first: a long shape must not overflow the cell index
	const double ka = L.origin_k + std::ceil((za - L.off_z * L.dz) / L.dz - EPS);
	const double kb = L.origin_k + std::floor((zb - L.off_z * L.dz) / L.dz + EPS) + 1;
	const int k0 = int(std::max(ka, double(L.k0)));
	const int k1 = int(std::min(kb, double(L.k1)));
	if (k0 < k1)
		row.push_back({k0, k1});
}

Shape make_shape(Shape::Node n) { return Shape(std::make_shared<const Shape::Node>(std::move(n))); }

std::pair<double, double> rotate(std::pair<double, double> p, double c, double s, double cx, double cz)
{
	const double x = p.first - cx, z = p.second - cz;
	return {cx + c * x - s * z, cz + s * x + c * z};
}
} // namespace

Row Shape::row(double x, const Lattice& L) const
{
	Row out;
	if (!_node)
		return out;
	const Node& n = *_node;
	switch (n.kind)
	{
	case Node::Rect:
		if (x >= n.x0 - EPS * L.dx && x <= n.x1 + EPS * L.dx)
			add_interval(out, n.z0, n.z1, L);
		break;
	case Node::Circle:
	{
		const double d = x - n.cx, h2 = n.r * n.r - d * d;
		if (h2 >= -EPS * L.dx * n.r)
		{
			const double h = std::sqrt(std::max(h2, 0.0));
			add_interval(out, n.cz - h, n.cz + h, L);
		}
		break;
	}
	case Node::Polygon:
	{
		std::vector<double> zs;
		const size_t m = n.pts.size();
		for (size_t e = 0; e < m; ++e)
		{
			const auto& p = n.pts[e];
			const auto& q = n.pts[(e + 1) % m];
			if ((p.first <= x && x < q.first) || (q.first <= x && x < p.first))
				zs.push_back(p.second + (x - p.first) * (q.second - p.second) / (q.first - p.first));
		}
		std::sort(zs.begin(), zs.end());
		Row raw;
		for (size_t c = 0; c + 1 < zs.size(); c += 2)
			add_interval(raw, zs[c], zs[c + 1], L);
		for (const auto& r : raw) // intervals of one polygon may touch after rounding
			out = geom::unite(out, Row{r});
		break;
	}
	case Node::Union:
		return geom::unite(n.a.row(x, L), n.b.row(x, L));
	case Node::Difference:
		return geom::subtract(n.a.row(x, L), n.b.row(x, L));
	case Node::Intersection:
		return geom::intersect(n.a.row(x, L), n.b.row(x, L));
	}
	return out;
}

Layer Shape::render(const Lattice& L) const
{
	Layer layer;
	layer.i0 = L.i0;
	layer.rows.reserve(std::max(0, L.i1 - L.i0));
	for (int i = L.i0; i < L.i1; ++i)
		layer.rows.push_back(row(L.x(i), L));
	return layer;
}

Shape Shape::rotated(double theta, double cx, double cz) const
{
	if (!_node || theta == 0)
		return *this;
	const Node& n = *_node;
	const double c = std::cos(theta), s = std::sin(theta);
	Node m = n;
	switch (n.kind)
	{
	case Node::Rect:
		m.kind = Node::Polygon;
		m.pts = {{n.x0, n.z0}, {n.x1, n.z0}, {n.x1, n.z1}, {n.x0, n.z1}};
		for (auto& p : m.pts)
			p = rotate(p, c, s, cx, cz);
		break;
	case Node::Circle:
		std::tie(m.cx, m.cz) = rotate({n.cx, n.cz}, c, s, cx, cz);
		break;
	case Node::Polygon:
		for (auto& p : m.pts)
			p = rotate(p, c, s, cx, cz);
		break;
	default:
		m.a = n.a.rotated(theta, cx, cz);
		m.b = n.b.rotated(theta, cx, cz);
	}
	return make_shape(m);
}

Shape Shape::translated(double dx, double dz) const
{
	if (!_node || (dx == 0 && dz == 0))
		return *this;
	Node m = *_node;
	switch (m.kind)
	{
	case Node::Rect:
		m.x0 += dx, m.x1 += dx, m.z0 += dz, m.z1 += dz;
		break;
	case Node::Circle:
		m.cx += dx, m.cz += dz;
		break;
	case Node::Polygon:
		for (auto& p : m.pts)
			p.first += dx, p.second += dz;
		break;
	default:
		m.a = m.a.translated(dx, dz);
		m.b = m.b.translated(dx, dz);
	}
	return make_shape(m);
}

Shape Shape::mirrored_x(double cx) const
{
	if (!_node)
		return *this;
	Node m = *_node;
	switch (m.kind)
	{
	case Node::Rect:
	{
		const double x0 = m.x0;
		m.x0 = 2 * cx - m.x1, m.x1 = 2 * cx - x0;
		break;
	}
	case Node::Circle:
		m.cx = 2 * cx - m.cx;
		break;
	case Node::Polygon:
		for (auto& p : m.pts)
			p.first = 2 * cx - p.first;
		break;
	default:
		m.a = m.a.mirrored_x(cx);
		m.b = m.b.mirrored_x(cx);
	}
	return make_shape(m);
}

double Shape::area() const
{
	if (!_node)
		return 0;
	const Node& n = *_node;
	switch (n.kind)
	{
	case Node::Rect:
		return (n.x1 - n.x0) * (n.z1 - n.z0);
	case Node::Circle:
		return M_PI * n.r * n.r;
	case Node::Polygon:
	{
		double a = 0;
		for (size_t e = 0; e < n.pts.size(); ++e)
		{
			const auto& p = n.pts[e];
			const auto& q = n.pts[(e + 1) % n.pts.size()];
			a += p.first * q.second - q.first * p.second;
		}
		return std::fabs(0.5 * a);
	}
	default:
		return -1;
	}
}

Shape box(double x0, double z0, double x1, double z1)
{
	Shape::Node n;
	n.kind = Shape::Node::Rect;
	if (x0 > x1 || z0 > z1)
		return Shape(); // inverted bounds: empty, as a box test lo <= v <= hi would be
	n.x0 = x0, n.x1 = x1, n.z0 = z0, n.z1 = z1;
	return make_shape(n);
}

Shape rect(double cx, double cz, double width_x, double width_z)
{
	return box(cx - 0.5 * width_x, cz - 0.5 * width_z, cx + 0.5 * width_x, cz + 0.5 * width_z);
}

Shape circle(double cx, double cz, double radius)
{
	Shape::Node n;
	n.kind = Shape::Node::Circle;
	n.cx = cx, n.cz = cz, n.r = radius;
	return make_shape(n);
}

Shape polygon(std::vector<std::pair<double, double>> xz)
{
	Shape::Node n;
	n.kind = Shape::Node::Polygon;
	n.pts = std::move(xz);
	return make_shape(n);
}

Shape rounded_rect(double cx, double cz, double width_x, double width_z, double r)
{
	if (r <= 0)
		return rect(cx, cz, width_x, width_z);
	const double hx = 0.5 * width_x - r, hz = 0.5 * width_z - r;
	return rect(cx, cz, width_x, 2 * hz) + rect(cx, cz, 2 * hx, width_z) + circle(cx - hx, cz - hz, r) +
		   circle(cx + hx, cz - hz, r) + circle(cx - hx, cz + hz, r) + circle(cx + hx, cz + hz, r);
}

Shape rounded_polygon(const std::vector<std::pair<double, double>>& p, double r)
{
	if (r <= 0)
		return polygon(p);
	const size_t n = p.size();
	std::vector<std::pair<double, double>> cut; // the two tangent points of each corner
	Shape discs;
	for (size_t i = 0; i < n; ++i)
	{
		const auto [vx, vz] = p[i];
		auto unit = [&](std::pair<double, double> q) {
			const double ux = q.first - vx, uz = q.second - vz, l = std::hypot(ux, uz);
			return std::pair{ux / l, uz / l};
		};
		const auto [ax, az] = unit(p[(i + n - 1) % n]);
		const auto [bx, bz] = unit(p[(i + 1) % n]);
		const double half = 0.5 * std::acos(std::clamp(ax * bx + az * bz, -1.0, 1.0));
		const double t = r / std::tan(half), c = r / std::sin(half), bl = std::hypot(ax + bx, az + bz);
		cut.push_back({vx + t * ax, vz + t * az});
		cut.push_back({vx + t * bx, vz + t * bz});
		discs = discs + circle(vx + c * (ax + bx) / bl, vz + c * (az + bz) / bl, r);
	}
	return polygon(cut) + discs;
}

Shape inflated_polygon(const std::vector<std::pair<double, double>>& p, double r)
{
	if (r <= 0)
		return polygon(p);
	const size_t n = p.size();
	double twice_area = 0; // > 0 for counter-clockwise order
	for (size_t i = 0; i < n; ++i)
		twice_area += p[i].first * p[(i + 1) % n].second - p[(i + 1) % n].first * p[i].second;
	const double out = twice_area > 0 ? 1.0 : -1.0;
	auto normal = [&](size_t e) { // outward unit normal of edge e = (p[e], p[e+1]), times r
		const double ex = p[(e + 1) % n].first - p[e].first, ez = p[(e + 1) % n].second - p[e].second, l = std::hypot(ex, ez);
		return std::pair{out * r * ez / l, -out * r * ex / l};
	};
	std::vector<std::pair<double, double>> hull; // each vertex moved out along both of its edges' normals
	Shape discs;
	for (size_t i = 0; i < n; ++i)
	{
		const auto [ax, az] = normal((i + n - 1) % n);
		const auto [bx, bz] = normal(i);
		hull.push_back({p[i].first + ax, p[i].second + az});
		hull.push_back({p[i].first + bx, p[i].second + bz});
		discs = discs + circle(p[i].first, p[i].second, r);
	}
	return polygon(hull) + discs;
}

namespace
{
Shape combine_shapes(Shape::Node::Kind k, const Shape& a, const Shape& b)
{
	Shape::Node n;
	n.kind = k;
	n.a = a, n.b = b;
	return make_shape(n);
}
} // namespace

Shape operator+(const Shape& a, const Shape& b)
{
	if (a.empty())
		return b;
	if (b.empty())
		return a;
	return combine_shapes(Shape::Node::Union, a, b);
}

Shape operator-(const Shape& a, const Shape& b)
{
	if (a.empty() || b.empty())
		return a;
	return combine_shapes(Shape::Node::Difference, a, b);
}

Shape operator&(const Shape& a, const Shape& b)
{
	if (a.empty() || b.empty())
		return Shape();
	return combine_shapes(Shape::Node::Intersection, a, b);
}

// solids

struct Solid::Node
{
	enum Kind { Extrusion, Sphere, CylinderZ, CylinderX, Sliced, Moved, Union, Difference, Intersection } kind = Extrusion;
	std::function<Shape(double)> at;
	Shape shape;
	double y0 = 0, y1 = 0;                // extrusion; cylinders: axial extent
	double cx = 0, cy = 0, cz = 0, r = 0; // sphere, cylinders; moved: shift (cx, cy, cz)
	double theta = 0, ox = 0, oz = 0;     // moved: then turned by theta about (ox, oz)
	Solid a, b;
};

namespace
{
Solid make_solid(Solid::Node n) { return Solid(std::make_shared<const Solid::Node>(std::move(n))); }
}

Shape Solid::slice(double y) const
{
	if (!_node)
		return Shape();
	const Node& n = *_node;
	switch (n.kind)
	{
	case Node::Extrusion:
		return (y >= n.y0 - EPS * std::fabs(n.y1 - n.y0) && y <= n.y1 + EPS * std::fabs(n.y1 - n.y0)) ? n.shape : Shape();
	case Node::Sphere:
	{
		const double h2 = n.r * n.r - (y - n.cy) * (y - n.cy);
		return h2 >= 0 ? circle(n.cx, n.cz, std::sqrt(h2)) : Shape();
	}
	case Node::CylinderZ:
	{
		const double h2 = n.r * n.r - (y - n.cy) * (y - n.cy);
		return h2 >= 0 ? box(n.cx - std::sqrt(h2), n.y0, n.cx + std::sqrt(h2), n.y1) : Shape();
	}
	case Node::CylinderX:
	{
		const double h2 = n.r * n.r - (y - n.cy) * (y - n.cy);
		return h2 >= 0 ? box(n.y0, n.cz - std::sqrt(h2), n.y1, n.cz + std::sqrt(h2)) : Shape();
	}
	case Node::Sliced:
		return n.at(y);
	case Node::Moved:
		return n.a.slice(y - n.cy).translated(n.cx, n.cz).rotated(n.theta, n.ox, n.oz);
	case Node::Union:
		return n.a.slice(y) + n.b.slice(y);
	case Node::Difference:
		return n.a.slice(y) - n.b.slice(y);
	case Node::Intersection:
		return n.a.slice(y) & n.b.slice(y);
	}
	return Shape();
}

Solid extrude(const Shape& s, double y0, double y1)
{
	if (y0 > y1)
		return Solid();
	Solid::Node n;
	n.kind = Solid::Node::Extrusion;
	n.shape = s, n.y0 = y0, n.y1 = y1;
	return make_solid(n);
}

Solid sphere(double cx, double cy, double cz, double radius)
{
	if (radius <= 0)
		return Solid();
	Solid::Node n;
	n.kind = Solid::Node::Sphere;
	n.cx = cx, n.cy = cy, n.cz = cz, n.r = radius;
	return make_solid(n);
}

Solid cylinder_z(double cx, double cy, double radius, double z0, double z1)
{
	if (radius <= 0)
		return Solid();
	Solid::Node n;
	n.kind = Solid::Node::CylinderZ;
	n.cx = cx, n.cy = cy, n.r = radius, n.y0 = z0, n.y1 = z1;
	return make_solid(n);
}

Solid cylinder_x(double cy, double cz, double radius, double x0, double x1)
{
	if (radius <= 0)
		return Solid();
	Solid::Node n;
	n.kind = Solid::Node::CylinderX;
	n.cy = cy, n.cz = cz, n.r = radius, n.y0 = x0, n.y1 = x1;
	return make_solid(n);
}

Solid Solid::translated(double dx, double dy, double dz) const
{
	if (!_node)
		return *this;
	Node n;
	n.kind = Node::Moved;
	n.a = *this, n.cx = dx, n.cy = dy, n.cz = dz;
	return make_solid(n);
}

Solid Solid::rotated(double theta, double cx, double cz) const
{
	if (!_node || theta == 0)
		return *this;
	Node n;
	n.kind = Node::Moved;
	n.a = *this, n.theta = theta, n.ox = cx, n.oz = cz;
	return make_solid(n);
}

Solid sliced(std::function<Shape(double y)> at)
{
	Solid::Node n;
	n.kind = Solid::Node::Sliced;
	n.at = std::move(at);
	return make_solid(n);
}

namespace
{
Solid combine_solid(Solid::Node::Kind k, const Solid& a, const Solid& b)
{
	Solid::Node n;
	n.kind = k;
	n.a = a, n.b = b;
	return make_solid(n);
}
} // namespace

Solid operator+(const Solid& a, const Solid& b) { return a.empty() ? b : b.empty() ? a : combine_solid(Solid::Node::Union, a, b); }
Solid operator-(const Solid& a, const Solid& b) { return (a.empty() || b.empty()) ? a : combine_solid(Solid::Node::Difference, a, b); }
Solid operator&(const Solid& a, const Solid& b) { return (a.empty() || b.empty()) ? Solid() : combine_solid(Solid::Node::Intersection, a, b); }

} // namespace geom
