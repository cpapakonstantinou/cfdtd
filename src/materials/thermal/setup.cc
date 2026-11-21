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
		fdtd::io::log("thermal", "loaded %s", strS.c_str());
		for (timestep = 1; timestep <= nr_max_time_steps_pump; timestep++) {
			for (j = 0; j <= subdomain.cells_y; j++) {
				temp_save[0][timestep][j] = double(atof(fdtd::io::read_text_entry(pfin).c_str()));
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
	MPI_Barrier(MPI_COMM_WORLD);

	if (fdtd::grid::world_rank == ROOTPROCESS) {
		fdtd::io::log("thermal", "temperature field loaded");
	}

	for (i = 0; i <= subdomain.inner_dim_x; ++i) {
		for (j = 0; j <= subdomain.inner_dim_y; ++j) {
			for (k = 0; k <= subdomain.inner_dim_z; ++k) {
				electron_temp[i + subdomain.inner_org_x][j + subdomain.inner_org_y][k + subdomain.inner_org_z] = electron_temp_tf[i][j][k];
			}
		}
	}
}

