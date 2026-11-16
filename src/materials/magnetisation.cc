#include "grid/subdomain.h"
#include "materials/state.h"
#include "grid/state.h"
#include <string>
#include "io/text.h"
#include "materials/materials.h"
#include "grid/lattice.h"
using fdtd::grid::subdomain;

namespace
{
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace

void fdtd::materials::read_magnetisation_field(const std::string &input_dir)
{
	int i;
	int j;
	int k;
	double valx;
	double valy;
	double valz;
	std::string strS = input_dir + "/Magnetization_x.txt"; //files each must have Nx+2 x Ny+2 x Nz/2 values
	std::fstream pfin;
	pfin.open(strS.c_str(), std::ios::in);
	std::string strS1 = input_dir + "/Magnetization_y.txt";
	std::fstream pfin1;
	pfin1.open(strS1.c_str(), std::ios::in);
	std::string strS2 = input_dir + "/Magnetization_z.txt";
	std::fstream pfin2;
	pfin2.open(strS2.c_str(), std::ios::in);
	if (pfin.is_open())
	{
		for (i = 0; i <= lat.cells[0] + 1; ++i) //extra cell required for interpolation
		{
			for (j = 0; j <= lat.cells[1] + 1; ++j)
			{
				for (k = 0; k <= lat.cells[2] + 1; ++k)
				{
					valx = double(atof(fdtd::io::read_text_entry(pfin).c_str()));
					valy = double(atof(fdtd::io::read_text_entry(pfin1).c_str()));
					valz = double(atof(fdtd::io::read_text_entry(pfin2).c_str()));

					if(i >= subdomain.start_x && i <= subdomain.start_x + subdomain.cells_x + 1
						&& j >= subdomain.start_y && j <= subdomain.start_y + subdomain.cells_y + 1
						&& k >= subdomain.start_z && k <= subdomain.start_z + subdomain.cells_z + 1){
							magnetization_x[i - subdomain.start_x][j - subdomain.start_y][k - subdomain.start_z] = valx;
							magnetization_y[i - subdomain.start_x][j - subdomain.start_y][k - subdomain.start_z] = valy;
							magnetization_z[i - subdomain.start_x][j - subdomain.start_y][k - subdomain.start_z] = valz;
					}
				}
			}
		}
		pfin.close();
		pfin1.close();
		pfin2.close();
	}
	return;
}
