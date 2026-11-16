#include "grid/subdomain.h"
#include "grid/state.h"
#include "materials/thermal/state.h"
#include "materials/thermal/thermal.h"
#include "io/text.h"
#include "materials/spec.h"
#include "grid/lattice.h"
using fdtd::constants::DEFAULT_TEMPERATURE;
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;

namespace
{
char _cstring[100]; ///< scratch buffer for a file name
} // namespace
using namespace fdtd::materials::thermal; // transitional: module state, thermal/state.h

void fdtd::materials::thermal::load_probe_temperature_1d(){

	int j;
	int timestep;
	sprintf(_cstring, "../OUTPUTPumpTemp_%d.txt", fdtd::grid::world_rank);
	std::string strS = _cstring;
	std::fstream pfin;
	pfin.open(strS.c_str(), std::ios::in);
	if (pfin.is_open()) {
		fdtd::io::log("thermal", "%s\t%s\n", strS.c_str(), "LOADED");
		for (timestep = 1; timestep <= nr_max_time_steps_pump; timestep++) {
			for (j = 0; j <= subdomain.cells_y; j++) {
				temp_save[0][timestep][j] = double(atof(fdtd::io::read_text_entry(pfin).c_str()));
				if (fdtd::grid::world_rank == 30 && timestep == 100000)
					fdtd::io::log("thermal", "Temp Here = %e\n", temp_save[0][timestep][j]);
			}
		}
	}
	else {
		for (timestep = 1; timestep <= nr_max_time_steps_pump; timestep++) {
			for (j = 0; j <= subdomain.cells_y; j++) {
				temp_save[0][timestep][j] = DEFAULT_TEMPERATURE;
			}
		}
	}

}

void fdtd::materials::thermal::load_static_temperature(){
	int i;
	int j;
	int k;
	MPI_Barrier(MPI_COMM_WORLD);
	fdtd::materials::thermal::read_temperature();
	MPI_Barrier(MPI_COMM_WORLD);
	fdtd::io::log("thermal", "LOADING: %d\t%d\t%d\n", fdtd::grid::lattice().inner(0), fdtd::grid::lattice().inner(1), fdtd::grid::lattice().inner(2));
	MPI_Barrier(MPI_COMM_WORLD);

	if (fdtd::grid::world_rank == ROOTPROCESS) {
		fdtd::io::log("thermal", "Temperature Field Loaded!\n");
	}

	for (i = 0; i <= subdomain.inner_dim_x; ++i) {
		for (j = 0; j <= subdomain.inner_dim_y; ++j) {
			for (k = 0; k <= subdomain.inner_dim_z; ++k) {
				electron_temp[i + subdomain.inner_org_x][j + subdomain.inner_org_y][k + subdomain.inner_org_z] = electron_temp_tf[i][j][k];
			}
		}
	}
}

