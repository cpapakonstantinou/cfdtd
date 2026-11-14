/**
 * @file app/structures.h
 * @brief Draw the structure of the run, background, blocks, the [geometry] shape, into the material maps
 */
#ifndef FDTD_APP_STRUCTURES_H
#define FDTD_APP_STRUCTURES_H

#include "geometry/anapole.h"
#include "geometry/block.h"

namespace fdtd::app
{
extern geom::AnapoleSpec anapole_input;            ///< [anapole] (apply_sim_config)
extern std::vector<geometry::Block> blocks;        ///< [block.*] in file order (apply_sim_config)
void create_structures();
} // namespace fdtd::app

#endif
