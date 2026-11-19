/**
 * @file materials/material_map.h
 * @brief The material maps of this rank: the code of every Yee E point, and cell centre
 *
 * A map holds a material slot per point, materials/slots.h. The geometry writes the inner maps,
 * the update scans and the monitors read them. The full maps cover every cell of the rank, PML
 * included; the mode source and its CPML read them.
 * @verbatim
      x, y, z            Ex, Ey, Ez points of the rank's inner block   [u][w][k], u = i - inner_org_x
      centre             cell centres, rank cells + 2, hydrodynamic,   [i][j][k]
      full_x, _y, _z     Ex, Ey, Ez points of every rank cell           [i][j][k]
   @endverbatim
 */
#ifndef FDTD_MATERIALS_MATERIAL_MAP_H
#define FDTD_MATERIALS_MATERIAL_MAP_H

#include "core/array.h"
#include "materials/material_field.h"
#include <vector>

namespace fdtd::materials
{
class MaterialMap
{
public:
	array3<int> x; ///< inner block, per E component
	array3<int> y; ///< inner block, per E component
	array3<int> z; ///< inner block, per E component
	array3<int> centre;                 ///< cell centres
	array3<int> full_x; ///< every rank cell, per E component
	array3<int> full_y; ///< every rank cell, per E component
	array3<int> full_z; ///< every rank cell, per E component

	/// the inner maps and the cell-centre map of the rank, inner block size, rank cell counts
	void allocate_inner(int ni, int nj, int nk, int cells_x, int cells_y, int cells_z);
	/// the full maps, rank cell counts
	void allocate_full(int cells_x, int cells_y, int cells_z);
	/// the full maps from the inner ones: a cell outside the inner block takes the slot of the
	/// nearest inner cell, so a structure that reaches the domain edge continues through the PML;
	/// `background` is the slot of a rank with no inner cells, per y index of the rank
	void extend_into_pml(int org_x, int org_y, int org_z, int ni, int nj, int nk, bool no_inner, const std::vector<int> &background);
	/// bounding rows (i, j) of every dispersive slot on the inner maps, in rank cell indices;
	/// i1 < 0: not on this rank
	std::vector<CellBox> dispersive_boxes(int org_x, int org_y) const;
	void release();
};

/// this rank's material maps
MaterialMap &material_map();

/// the maps by their solver names
inline array3<int> &matind_x = material_map().x;
inline array3<int> &matind_y = material_map().y;
inline array3<int> &matind_z = material_map().z;
inline array3<int> &matind = material_map().centre;
inline array3<int> &matind_x_full = material_map().full_x;
inline array3<int> &matind_y_full = material_map().full_y;
inline array3<int> &matind_z_full = material_map().full_z;
} // namespace fdtd::materials

#endif
