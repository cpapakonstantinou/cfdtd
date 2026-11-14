/**
 * @file app/release.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @date 2022-06-22
 *
 * @brief Contains a single function fdtd::app::release_run_memory() which releases the memory for various fields and monitors, 
 * which have been allocated by functions called in the initialization of the simulation before
 * see memory_allocation.cpp. It will be called at the end of each simulation.
 *
 * @copyright Copyright 2022
 */
#include "materials/state.h"
#include "materials/material_map.h"
#include "app/run_config.h"
#include "monitors/state.h"
#include "materials/thermal/state.h"
#include "sources/state.h"
#include "sources/aux_grid.h"
#include "materials/hydro/state.h"
#include "app/setup.h"
#include "monitors/volume_dft.h"
#include "monitors/field_trackers.h"
#include "monitors/flux_boxes.h"
#include "monitors/art_planes.h"
#include "fields/fields.h"
#include "boundary/cpml.h"
#include "io/binary.h"
using fdtd::fields::ex_nl;
using fdtd::fields::ey_nl;
using fdtd::fields::ez_nl;
using fdtd::fields::hx_nl;
using fdtd::fields::hy_nl;
using fdtd::fields::hz_nl;
using fdtd::materials::DrudeCP;
using fdtd::materials::Lorentz;
using fdtd::materials::matind_x;
using fdtd::materials::matind_y;
using fdtd::materials::matind_z;
using fdtd::materials::nr_disp_mat_d;

namespace
{
const fdtd::app::RunConfig &run_cfg = fdtd::app::run_config(); ///< the run configuration
} // namespace

/**
 * @brief Releases all the memory, of all fields and monitors, which have been used in the simulation
 */
void fdtd::app::release_run_memory()
{
	fdtd::app::material_state_intact();

	fdtd::fields::fields().release(); // E, H and their history
	if (fdtd::materials::thermal::temp_save != NULL)
		fdtd::materials::thermal::temp_save.reset();
	if (fdtd::materials::thermal::electron_temp != NULL)
		fdtd::materials::thermal::electron_temp.reset();
	if (fdtd::materials::thermal::lattice_temp != NULL)
		fdtd::materials::thermal::lattice_temp.reset();
	if (fdtd::monitors::energy_density != NULL)
		fdtd::monitors::energy_density.reset();

	fdtd::materials::thermal::electron_temp_tf.reset();
	fdtd::materials::hydro::nd_z_prev.reset();
	fdtd::materials::magnetization_x.reset(), fdtd::materials::magnetization_y.reset(), fdtd::materials::magnetization_z.reset();

	// the monitors' DFT arrays
	fdtd::monitors::field_trackers().release();
	fdtd::monitors::art_planes().release();
	fdtd::monitors::flux_boxes().release();
	fdtd::monitors::dft_planes().release();
	fdtd::monitors::volume_dft().release();

	fdtd::materials::dispersive_state().release(); // P, J, poles and the nonlocal / hydrodynamic arrays

	fdtd::materials::material_map().release();

	fdtd::boundary::cpml().release(); // CPML profiles and psi arrays

	fdtd::sources::aux_grid().release();

	return;
}