/**
 * @file geometry/render.cc
 * @brief Render geometry-layer drawings into the material maps; the plasmonic anapole
 */
#include "grid/subdomain.h"
#include "core/carray.h"
#include "geometry/csg.h"
#include "materials/component.h"
#include "materials/material_map.h"
#include "materials/codes.h"
#include "core/constants.h"
#include "io/text.h"
#include "geometry/render.h"
#include "grid/state.h"
#include "materials/nonlinear/state.h"
#include "grid/lattice.h"
#include "geometry/scene.h"
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;
using fdtd::materials::matind_x;
using fdtd::materials::matind_y;
using fdtd::materials::matind_z;
using fdtd::materials::nl_code;

namespace
{
const fdtd::geometry::Scene &sc = fdtd::geometry::scene();            ///< what the builders draw into
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace

void render_parts(const std::vector<geom::Part>& parts, int ci, int cj, int ck, Staircase staircase)
{
	const int si = int(subdomain.output_start[0]), sj = int(subdomain.output_start[1]), sk = int(subdomain.output_start[2]);
	const int ni = int(subdomain.output_size[0]), nj = int(subdomain.output_size[1]), nk = int(subdomain.output_size[2]);
	// lattice of one component: sample offsets (x, z) in cells
	auto lattice = [&](double ox, double oz) {
		geom::Lattice L;
		L.dx = lat.d[0], L.dz = lat.d[2], L.off_x = ox, L.off_z = oz, L.origin_i = ci, L.origin_k = ck;
		L.i0 = si, L.i1 = si + ni, L.k0 = sk, L.k1 = sk + nk;
		return L;
	};
	const geom::Lattice Lx = lattice(0.5, 0), Ly = lattice(0, 0), Lz = lattice(0, 0.5), Lc = lattice(0.5, 0.5);
	auto y_at = [&](int j, double off) { return (j - cj) * lat.d[1] + off * lat.d[1]; };
	auto fill = [&](ctensor<int>& map, int i, int j, const geom::Row& row, int material) {
		for (const auto& run : row)
		{
			int* p = &map[i - si][j - sj][run.k0 - sk];
			std::fill(p, p + (run.k1 - run.k0), material);
		}
	};

	for (const auto& part : parts)
		for (int j = sj; j < sj + nj; ++j)
		{
			const geom::Shape at_y = part.body.slice(y_at(j, 0));
			const geom::Shape at_yf = part.body.slice(y_at(j, 0.5));
			if (at_y.empty() && at_yf.empty())
				continue;
			for (int i = si; i < si + ni; ++i)
			{
				if (!sc.staircase)
				{
					fill(matind_x, i, j, at_y.row(Lx.x(i), Lx), part.material);
					fill(matind_y, i, j, at_yf.row(Ly.x(i), Ly), part.material);
					fill(matind_z, i, j, at_y.row(Lz.x(i), Lz), part.material);
					continue;
				}
				const geom::Row row = staircase == Staircase::CellCentre
										  ? at_yf.row(Lc.x(i), Lc)
										  : geom::intersect(geom::intersect(at_y.row(Lx.x(i), Lx), at_yf.row(Ly.x(i), Ly)),
															   at_y.row(Lz.x(i), Lz));
				for (auto* map : {&matind_x, &matind_y, &matind_z})
					fill(*map, i, j, row, part.material);
			}
		}
}

void create_plasmonic_anapole(const geom::AnapoleSpec &a, int material, int ci, int cj, int ck)
{
	// drawn about the patch centre cell (ci, cj, ck), with the patch material, nonlinear offset included
	const auto parts = geom::anapole_parts(a, material, 0, 0, 0);
	if (fdtd::grid::world_rank == ROOTPROCESS)
	{
		const double s_neck = (a.neck_half_length / a.patch_edge - a.lobe_x0) / a.lobe_a;
		fdtd::io::log("geometry", "Plasmonic anapole: patch %.4g nm, thickness %.4g nm, corner radius %.4g nm, neck %.4g x %.4g nm, lobe tip %.4g nm, "
					"neck-face lobe half-height %.4g nm\n",
					a.patch_edge * 1e9, a.thickness * 1e9, a.corner_radius * 1e9, 2 * a.neck_half_length * 1e9, 2 * a.neck_half_height * 1e9,
					a.patch_edge * (a.lobe_x0 + a.lobe_a) * 1e9,
					a.patch_edge * a.lobe_b * std::pow(s_neck, a.lobe_p) * std::sqrt(1 - s_neck * s_neck) * 1e9);
		if (a.dimer_material >= 0 && (a.dimer_radius < lat.d[0] || a.dimer_gap < lat.d[0]))
			fdtd::io::log("geometry", "  warning: dimer features are below the grid step (%.3g nm)\n", lat.d[0] * 1e9);
	}
	render_parts(parts, ci, cj, ck, Staircase::CellCentre);
}
