/**
 * @file materials/thermal/output.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-06-15
 * 
 * @brief Houses parrallel read-in and write-out functions
 * 
 * @details Contains functions to write out temperature data for the
 * whole domain. As well, there is a function to read-in 3D Temperature data
 * in parallel
 * 
 * @copyright Copyright 2022
 * 
 */
#include "grid/state.h"
#include "materials/thermal/state.h"
#include "materials/thermal/thermal.h"
#include "grid/subdomain.h"
#include "io/h5io.h"

#include <vector>
#include "materials/spec.h"
#include "grid/lattice.h"
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;

namespace
{
} // namespace
using namespace fdtd::materials::thermal; // transitional: module state, thermal/state.h

namespace
{
/// This rank's, subdomain.output_size^3, block of the, globalsizes^3, temperature grid. Neighbouring
/// blocks share one boundary plane, as in the field arrays.
h5io::Block temperature_block()
{
	const hsize_t global[3] = {hsize_t((subdomain.cpml_only ? 0L : long(fdtd::grid::lattice().inner(0) + 1))), hsize_t((subdomain.cpml_only ? 0L : long(fdtd::grid::lattice().inner(1) + 1))), hsize_t((subdomain.cpml_only ? 0L : long(fdtd::grid::lattice().inner(2) + 1)))};
	const hsize_t offset[3] = {hsize_t(subdomain.output_start[0]), hsize_t(subdomain.output_start[1]), hsize_t(subdomain.output_start[2])};
	const hsize_t count[3] = {hsize_t(subdomain.output_size[0]), hsize_t(subdomain.output_size[1]), hsize_t(subdomain.output_size[2])};
	return h5io::Block::volume(global, offset, count);
}
} // namespace

/// Writes the electron temperature to <output_dir>/electron_temperature_<timestep>.h5
/// dataset /electron_temperature, attribute domain_size.
void fdtd::materials::thermal::save_temperature(int step)
{
	if (subdomain.cpml_only)
		return;
	h5io::File file(str_outfilepathname_t + std::to_string(step) + str_fileextension, fdtd::grid::output_comm);
	file.attribute("domain_size", std::vector<int>{fdtd::grid::lattice().inner(0), fdtd::grid::lattice().inner(1), fdtd::grid::lattice().inner(2)});
	file.attribute("timestep", std::vector<int>{step});
	h5io::write(file, "electron_temperature", temperature_block(), &electron_temp_tf[0][0][0]);
}

/// Loads the initial electron temperature from <input_dir>/electron_temperature.h5,
/// a file written by SaveOutputTemperature on the same domain.
void fdtd::materials::thermal::read_temperature()
{
	if (subdomain.cpml_only)
		return;
	h5io::File file(str_infilepathname_t + str_fileextension, fdtd::grid::output_comm, h5io::File::Access::read);
	const std::vector<int> size = file.int_attribute("domain_size");
	if (size != std::vector<int>{fdtd::grid::lattice().inner(0), fdtd::grid::lattice().inner(1), fdtd::grid::lattice().inner(2)})
	{
		if (fdtd::grid::world_rank == ROOTPROCESS)
			fprintf(stderr, "error: %s: domain_size does not match this simulation\n", file.path().c_str());
		MPI_Abort(MPI_COMM_WORLD, 4);
	}
	h5io::read(file, "electron_temperature", temperature_block(), &electron_temp_tf[0][0][0]);
}
