/**
 * @file geometry/scene.h
 * @brief What the builders need to know about the run besides their own dimensions
 *
 * Set by the application before it draws the structure, app::create_structures; the builders read
 * it instead of the run configuration, docs/dev/Roadmap.md §2.5.
 */
#ifndef FDTD_GEOMETRY_SCENE_H
#define FDTD_GEOMETRY_SCENE_H

namespace fdtd::geometry
{
struct Scene
{
	int shape = 0;              ///< the [geometry] shape, shape codes, geometry/shapes.h; a few builders draw variants by it
	bool draw = false;          ///< a shape is configured, shape != none; the antenna builders draw nothing otherwise
	bool staircase = false;     ///< one material for the three components of a cell, [grid] staircase
	int medium[3] = {};         ///< background media below, in and above the middle layer, material codes
	int centre[3] = {};         ///< the structure's centre cell
	int gap_material = 0;       ///< the gap fill of round_hybrid, material code, nonlinear offset included
};
const Scene &scene();
void set_scene(const Scene &s);
} // namespace fdtd::geometry

#endif
