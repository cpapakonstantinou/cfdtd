#pragma once
/**
 * @file geometry/csg.h
 * @brief Draw-then-render geometry: shapes are combined with set operators, then rendered once.
 *
 * Planes and axes: drawings live in the x-z plane, the lithography plane; y is the
 * extrusion axis. A material array is [i][j][k] with k (z) contiguous, so a raster row is
 * one x index i holding sorted runs of k.
 *
 * Rendering rule: a cell takes the material when its sample point is inside the shape
 * boundary included. The sample point of a Yee component is offset by half a cell along
 * that component's axis. Point sampling commutes with union, difference and intersection, so
 * the set operations run on the integer runs of each primitive, in O(runs) per row, and the
 * result is exactly the sampling of the combined shape.
 *
 * 3D shapes are stacks of 2D drawings along y: an extrusion repeats one drawing over a y
 * range, a sphere is a circle per slice.
 */
#include <functional>
#include <memory>
#include <vector>

namespace geom
{

/// one row of a raster: sorted, disjoint half-open runs [k0, k1)
struct Run
{
	int k0, k1;
	bool operator==(const Run& o) const { return k0 == o.k0 && k1 == o.k1; }
};
using Row = std::vector<Run>;

Row unite(const Row& a, const Row& b);
Row subtract(const Row& a, const Row& b);
Row intersect(const Row& a, const Row& b);

/**
 * @brief Where the samples of one plane sit
 * @verbatim
     z ^    k+1 .     .     .
       |     k  .     x     .        x = sample (i, k) of a component with
       |   k-1  .     .     .            off_x = 1/2, off_z = 0 (Ex)
       |       i-1    i    i+1
       +---------------------> x
   @endverbatim
 */
/// sample lattice of one plane: sample (i, k) sits at
/// x = (i - origin_i) dx + off_x dx,  z = (k - origin_k) dz + off_z dz
/// drawings are in metres relative to the origin cell; offsets are in cells
struct Lattice
{
	double dx = 1, dz = 1;
	double off_x = 0, off_z = 0;
	int origin_i = 0, origin_k = 0;
	int i0 = 0, i1 = 0; ///< rows rendered: [i0, i1)
	int k0 = 0, k1 = 0; ///< runs are clipped to [k0, k1)
	double x(int i) const { return (i - origin_i) * dx + off_x * dx; }
	double z(int k) const { return (k - origin_k) * dz + off_z * dz; }
};

/// raster of a 2D drawing on a lattice: rows[i - i0]
struct Layer
{
	int i0 = 0;
	std::vector<Row> rows;
	long cells() const;
};

/// a 2D drawing in the x-z plane [m]; build with the primitives below and the operators + - &
class Shape
{
public:
	struct Node;
	Shape() = default;
	explicit Shape(std::shared_ptr<const Node> n) : _node(std::move(n)) {}

	Layer render(const Lattice& lattice) const;
	Row row(double x, const Lattice& lattice) const; ///< the runs of one row, sampled at x
	Shape rotated(double theta, double cx, double cz) const; ///< rotate by theta about (cx, cz)
	Shape translated(double dx, double dz) const;          ///< move by (dx, dz)
	Shape mirrored_x(double cx) const;                     ///< mirror x -> 2 cx - x
	bool empty() const { return !_node; }
	/// exact area for a primitive; < 0 for a composite, use Layer::cells for the mesh area
	double area() const;

	const Node* node() const { return _node.get(); }

private:
	std::shared_ptr<const Node> _node;
};

/**
 * @brief Axis-aligned rectangle given by its centre and widths, empty for a negative width
 * @verbatim
     z ^
       |    +-----------+  ---
       |    |           |   ^
       |    |  (cx,cz)  |   | width_z
       |    |     +     |   v
       |    +-----------+  ---
       |    |<-width_x->|
       +---------------------> x
   @endverbatim
 */
Shape rect(double cx, double cz, double width_x, double width_z);

/**
 * @brief Axis-aligned rectangle given by its corners: [x0, x1] x [z0, z1]; empty when x0 > x1 or z0 > z1
 * @verbatim
     z ^
    z1 |    +-----------+
       |    |           |
    z0 |    +-----------+
       +----+-----------+----> x
            x0          x1
   @endverbatim
 */
Shape box(double x0, double z0, double x1, double z1);

/**
 * @brief Disc of radius r centred at (cx, cz)
 * @verbatim
     z ^     .-~~~-.
       |   /                |  |  (cx,cz)--|--> r
       |   \         /
       |     `-___-'
       +---------------> x
   @endverbatim
 */
Shape circle(double cx, double cz, double radius);

/**
 * @brief Closed polygon through the points (x, z), even-odd rule
 * @verbatim
     z ^        p2
       |        /\
       |       /  \          the last point connects back to the first
       |   p0 +----+ p1
       +---------------> x
   @endverbatim
 */
Shape polygon(std::vector<std::pair<double, double>> xz);

/**
 * @brief Rectangle with corners rounded by radius r: two rectangles plus four discs
 * @verbatim
     z ^     .---------.       ---
       |    /     r  .-|-.      ^
       |   |         | + |      |
       |   |    +    |   |      | width_z
       |   |         |   |      |
       |    \        `-.-'      v
       |     `---------'       ---
       |   |<-width_x->|
       +----------------------> x
   @endverbatim
 */
Shape rounded_rect(double cx, double cz, double width_x, double width_z, double r);

/**
 * @brief Convex polygon with every corner filleted by radius r
 * @verbatim
     z ^          .-.             each corner: a disc of radius r tangent to both edges,
       |         / c \            centre c at r / sin(phi/2) from the vertex along the
       |        /     \           bisector, phi = corner angle; the polygon through the
       |       /       \          tangent points t fills the rest
       |    t +.       .+ t
       |    .-'         `-.
       |   ( c   ------   c )
       |    `-+---------+-'
       +----------------------> x
   @endverbatim
 * r <= 0 gives the sharp polygon. r must fit the shortest edge, not checked.
 */
Shape rounded_polygon(const std::vector<std::pair<double, double>>& xz, double r);

/**
 * @brief Convex polygon grown by r in every direction: every point within r of it, Minkowski sum
 *        with a disc. Its bounding box is the polygon's plus r on each side.
 * @verbatim
     z ^     .-------------.          the edges move out by r along their normals;
       |    /   .-------.   \         each vertex becomes an arc of radius r
       |   |   /         \   |
       |   |  /  polygon  \  |
       |   | +-------------+ |
       |    \               /
       |     `-------------'
       +----------------------> x
   @endverbatim
 * rounded_polygon(p, r) is inflated_polygon of p moved in by r. r <= 0 gives the polygon.
 */
Shape inflated_polygon(const std::vector<std::pair<double, double>>& xz, double r);

Shape operator+(const Shape& a, const Shape& b); ///< union
Shape operator-(const Shape& a, const Shape& b); ///< difference
Shape operator&(const Shape& a, const Shape& b); ///< intersection

/// a 3D body: a stack of drawings along y [m]
class Solid
{
public:
	struct Node;
	Solid() = default;
	explicit Solid(std::shared_ptr<const Node> n) : _node(std::move(n)) {}
	/// the drawing of the slice at height y, empty Shape when the slice misses the body
	Shape slice(double y) const;
	Solid translated(double dx, double dy, double dz) const; ///< move by (dx, dy, dz)
	Solid rotated(double theta, double cx, double cz) const;  ///< turn by theta about the y axis through (cx, cz)
	bool empty() const { return !_node; }

private:
	std::shared_ptr<const Node> _node;
};

/**
 * @brief The drawing s, repeated for y0 <= y <= y1, empty when y0 > y1
 * @verbatim
     y ^
    y1 |   +=========+     every slice in [y0, y1] is the drawing s
       |   |    s    |
    y0 |   +=========+
       +-----------------> x, drawing in the x-z plane
   @endverbatim
 */
Solid extrude(const Shape& s, double y0, double y1);

/**
 * @brief Ball of radius r: the slice at height y is a disc of radius sqrt(r^2 - (y - cy)^2)
 * @verbatim
     y ^      .-~~-.        <- small disc
       |    /        \
    cy |   |    +     |     <- disc of radius r
       |    \        /
       |      `-__-'        <- small disc
       +-----------------> x
   @endverbatim
 */
Solid sphere(double cx, double cy, double cz, double radius); ///< empty for radius <= 0

/**
 * @brief Cylinder along z, axis at (cx, cy), from z0 to z1: each slice is the strip
 *        |x - cx| <= sqrt(r^2 - (y - cy)^2), z0 <= z <= z1, default: the whole domain. Empty for r <= 0.
 * @verbatim
     y ^   .-~~-.
       |  /      \        cross-section in the x-y plane
    cy | |   +----|-> r
       |  \      /
       |   `-__-'
       +--------------> x
   @endverbatim
 */
Solid cylinder_z(double cx, double cy, double radius, double z0 = -1e3, double z1 = 1e3);

/**
 * @brief Cylinder along x, axis at (cy, cz), from x0 to x1: each slice is the strip
 *        |z - cz| <= sqrt(r^2 - (y - cy)^2), x0 <= x <= x1. Empty for r <= 0.
 * @verbatim
     y ^   .-~~-.
       |  /      \        cross-section in the z-y plane
    cy | |   +----|-> r
       |  \      /
       |   `-__-'
       +--------------> z
   @endverbatim
 */
Solid cylinder_x(double cy, double cz, double radius, double x0, double x1);

/// A body given directly by its slices: `at(y)` returns the drawing at height y, empty when the
/// slice misses the body. For shapes whose outline varies with y, e.g. rounded edges.
Solid sliced(std::function<Shape(double y)> at);
Solid operator+(const Solid& a, const Solid& b);
Solid operator-(const Solid& a, const Solid& b);
Solid operator&(const Solid& a, const Solid& b);

} // namespace geom
