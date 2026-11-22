#pragma once
/**
 * @file geometry/render.h
 * @brief Render geometry-layer drawings into this rank's material maps
 */
#include "geometry/drawings.h"

/// How a staircased run, `[grid] staircase = true`, assigns one material to all three components
enum class Staircase
{
	AllComponents, ///< the cell takes the material when all three Yee sample points are inside
	CellCentre     ///< the cell takes the material when its centre is inside
};

/**
 * @brief Write each part's material into this rank's matind_x / _y / _z, in drawing order.
 *
 * The drawings are in metres relative to the origin cell (ci, cj, ck). A component takes the
 * material when its Yee sample point is inside:
 *
 *            Ex: (i + 1/2, j,       k      )
 *            Ey: (i,       j + 1/2, k      ), cells, relative to the origin cell
 *            Ez: (i,       j,       k + 1/2)
 *
 * When staircasing, one material is written to all three components according to `staircase`.
 */
void render_parts(const std::vector<geom::Part>& parts, int ci, int cj, int ck, Staircase staircase = Staircase::AllComponents);
