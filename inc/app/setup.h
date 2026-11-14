/**
 * @file app/setup.h
 * @brief The set-up and tear-down of one run, driver phases: allocations, initial state, release
 *
 * @verbatim
      setup_domain()                         grid sizes, CPML, cell size and time step from the config
      initialise_run()                       per-run constants, sources, material tables
      allocate_fields()                      E, H, history and model state arrays
      allocate_planes_and_material_maps()    animation flags, display planes, material maps
      allocate_aux_source()                  the 1D grid of the TF/SF source
      allocate_scattering_and_absorption()   flux boxes, volume DFT, ART planes
      release_run_memory()                   everything the run allocated, between sweeps and at the end
   @endverbatim
 * The allocators return false on failure.
 */
#ifndef FDTD_APP_SETUP_H
#define FDTD_APP_SETUP_H

namespace fdtd::app
{
void setup_domain();
bool initialise_run();
/// DT and 1/d from the configuration, no MPI; also used by `fdtd check` and `fdtd signal`
void set_time_step();
bool allocate_fields();
/// the dispersive state, after the geometry, rows of the materials' cells only where possible
bool allocate_material_state();
/// false, and a message, when the dispersive state was written outside a material box
bool material_state_intact();
bool allocate_planes_and_material_maps();
bool allocate_aux_source();
bool allocate_scattering_and_absorption();
void release_run_memory();
} // namespace fdtd::app

#endif
