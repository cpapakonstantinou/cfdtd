/**
 * @file materials/material_map.cc
 * @brief The material maps of this rank, materials/material_map.h
 */
#include "materials/material_map.h"
#include "materials/codes.h"
#include "core/views.h"
#include <algorithm>

namespace fdtd::materials
{
MaterialMap &material_map()
{
	static MaterialMap map;
	return map;
}

void MaterialMap::allocate_inner(int ni, int nj, int nk, int cells_x, int cells_y, int cells_z)
{
	x = field3<int>(ni, nj, nk);
	y = field3<int>(ni, nj, nk);
	z = field3<int>(ni, nj, nk);
	centre = field3<int>(cells_x + 2, cells_y + 2, cells_z + 2);
}

void MaterialMap::allocate_full(int cells_x, int cells_y, int cells_z)
{
	full_x = field3<int>(cells_x, cells_y, cells_z);
	full_y = field3<int>(cells_x, cells_y, cells_z);
	full_z = field3<int>(cells_x, cells_y, cells_z);
}

void MaterialMap::extend_into_pml(int org_x, int org_y, int org_z, int ni, int nj, int nk, bool no_inner, const std::vector<int> &background)
{
	ctensor<int> *inner[3] = {&x, &y, &z};
	ctensor<int> *full[3] = {&full_x, &full_y, &full_z};
	for (int c = 0; c < 3; ++c)
	{
		const view3<int> f = *full[c];
		const int cx = int(full[c]->extent(0));
		const int cy = int(full[c]->extent(1));
		const int cz = int(full[c]->extent(2));
		for (int i = 0; i < cx; ++i)
			for (int j = 0; j < cy; ++j)
				for (int k = 0; k < cz; ++k)
				{
					if (no_inner)
					{
						f[i][j][k] = background[std::clamp(j, 0, int(background.size()) - 1)];
						continue;
					}
					const view3<int> in = *inner[c];
					f[i][j][k] = in[std::clamp(i - org_x, 0, ni)][std::clamp(j - org_y, 0, nj)][std::clamp(k - org_z, 0, nk)];
				}
	}
}

std::vector<CellBox> MaterialMap::dispersive_boxes(int org_x, int org_y) const
{
	std::vector<CellBox> box(nr_disp_mat_d, CellBox{1 << 30, -1, 1 << 30, -1});
	for (const ctensor<int> *map : {&x, &y, &z})
	{
		const view3<int> v = *const_cast<ctensor<int> *>(map);
		if (!v)
			continue;
		for (size_t u = 0; u < map->extent(0); ++u)
			for (size_t w = 0; w < map->extent(1); ++w)
				for (size_t k = 0; k < map->extent(2); ++k)
				{
					int code = v[u][w][k];
					if (code >= nl_code)
						code -= nl_code;
					const int mm = code - DispMateOffsetD;
					if (mm < 0 || mm >= nr_disp_mat_d)
						continue;
					CellBox &b = box[mm];
					const int i = int(u) + org_x;
					const int j = int(w) + org_y;
					b.i0 = std::min(b.i0, i), b.i1 = std::max(b.i1, i), b.j0 = std::min(b.j0, j), b.j1 = std::max(b.j1, j);
				}
	}
	return box;
}

void MaterialMap::release()
{
	for (ctensor<int> *m : {&x, &y, &z, &centre, &full_x, &full_y, &full_z})
		m->reset();
}
} // namespace fdtd::materials
