/**
 * @file geometry/state.cc
 * @brief Definitions of the geometry module state, geometry/state.h
 */
#include "geometry/state.h"
#include "geometry/scene.h"

namespace fdtd::geometry
{
int j_offset_uni;
bool shg_bulk;
bool thg_bulk;
double corner_radius;
bool add_free_electron_layer;
double antenna_gap;
} // namespace fdtd::geometry

namespace fdtd::geometry
{
namespace
{
Scene the_scene; ///< set by the application before the structure is drawn
} // namespace
const Scene &scene() { return the_scene; }
void set_scene(const Scene &s) { the_scene = s; }
} // namespace fdtd::geometry
