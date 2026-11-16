
/**
 * @file grid/decomposition.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-01-19
 * 
 * @brief Ini file, housing initialization functions
 * 
 * @details This file is handles several things:\n
 * - The read-in of the varying ini file(s)
 * -  * and execution of the whole script. It walks through the creation of the
 * domains and arrays, and manages the FDTD loop with all proper updates to
 * appropriate monitors.
 * 
 * @warning Could perhaps be stricter to 'ini' functions
 * 
 * @copyright Copyright 2022
 */
#include <algorithm>
#include "grid/subdomain.h"
#include "grid/lattice.h"
#include <string>
#include "io/text.h"
#include "grid/state.h"
#include <cstdarg>
#include <vector>
#include <cstdio>
using fdtd::constants::ROOTPROCESS;

namespace
{
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::grid; // transitional: module state, grid/state.h

// every rank runs the decomposition; only the root reports

/**
 * @brief Performs the domain decomposition for parallel execution
 * @details It's exectuted only by Rootprocess, then the data is sent to all MPI processes
 * via MPI_broadcast
 * @param why set to the reason when a step fails
 *
 * @return true when an error has occured
 * @return false when the end of the function is reached successfully
 *
 */
/**
 * @brief This rank's CPML ranges along one axis, in local cell indices
 *
 * @verbatim
      global:  0 ........ nea [ physical domain ) far ........ n_tot
      rank:          [gmin ................. gmax)        n = gmax - gmin
      type 0, no CPML,          nea < gmin, gmax <= far
      type 1, only near CPML,   gmax <= nea
      type 2, only far CPML,    gmin > far
      type 3, near CPML + domain, type 4, domain + far CPML, type 5, both CPMLs
   @endverbatim
 * min / max are the local indices where the near CPML ends and the far CPML starts, nea - gmin,
 * far - gmin, outside [0, n] when the interface is not on this rank. A rank that is all CPML along
 * the axis has no inner cells, and its min / max describe an empty range beyond every local index
 * BELOW = -n or ABOVE = n + 100000: with the true indices, its cell on the interface, i = n or
 * 0, would pass the inner test, and that rank has no material maps, segfault in the E scans.
 * nea_cpml_end / far_cpml_start are the same indices as origins of the psi arrays and loop bounds: the
 * only-near rank's nea_cpml_end is n, all cells, the only-far rank's far_cpml_start is 0, psi from 0.
 * pad_nea / pad_far are the distances from the global start / end, the CPML grading index offset.
 * The nea_cpml_end / far_cpml_start values keep BELOW (-n) / ABOVE (n + 100000) where a face is not on this rank: they are
 * psi-array origins and loop bounds, not comparisons. Checked with rank edges exactly on the
 * interfaces, parallel.chunks_* = "0 8 44" of 52 cells, on each axis and all three together.
 */
struct AxisCpml
{
	int min;
	int max;
	int nea_cpml_end;
	int far_cpml_start;
	int pad_nea;
	int pad_far;
	int inner_dim;
	int inner_org;
	int type;
	int shift_output;
	bool only_cpml;
	bool no_cpml;
	void apply(int &mn, int &mx, int &mn_e, int &mx_e, int &p_n, int &p_f, int &dim, int &org, int &t, int &shift) const
	{
		mn = min, mx = max, mn_e = nea_cpml_end, mx_e = far_cpml_start, p_n = pad_nea, p_f = pad_far, dim = inner_dim, org = inner_org, t = type;
		if (!only_cpml)
			shift = shift_output;
	}
};

static AxisCpml axis_cpml_ranges(int gmin, int n, int nea, int far, int n_tot)
{
	const int gmax = gmin + n;
	const int BELOW = -n;
	const int ABOVE = n + 100000;
	const int a = nea - gmin;
	const int b = far - gmin;
	AxisCpml r{};
	r.min = a; // the interfaces in local indices, they may lie outside [0, n]
	r.max = b;
	r.pad_nea = gmin;
	r.pad_far = n_tot - gmax;
	r.only_cpml = gmax <= nea || gmin > far;
	if (gmin > far) // all far CPML: no inner cells, an empty inner range below every index
		r.min = BELOW, r.max = BELOW, r.nea_cpml_end = BELOW, r.far_cpml_start = 0, r.type = 2;
	else if (gmax <= nea) // all near CPML: an empty inner range above every index
		r.min = ABOVE, r.max = ABOVE, r.nea_cpml_end = n, r.far_cpml_start = ABOVE, r.type = 1;
	else if (gmin <= nea && gmax > far)
		r.nea_cpml_end = a, r.far_cpml_start = b, r.inner_dim = far - nea, r.inner_org = a, r.type = 5;
	else if (gmin <= nea)
		r.nea_cpml_end = a, r.far_cpml_start = ABOVE, r.inner_dim = n - a, r.inner_org = a, r.type = 3;
	else if (gmax <= far)
		r.nea_cpml_end = BELOW, r.far_cpml_start = ABOVE, r.inner_dim = n, r.inner_org = 0, r.type = 0, r.no_cpml = true;
	else
		r.nea_cpml_end = BELOW, r.far_cpml_start = b, r.inner_dim = b, r.inner_org = 0, r.type = 4;
	r.shift_output = gmin + r.inner_org - nea;
	return r;
}

namespace
{
/// The steps of set_lattice_data, in order. They share the rank's position, its extent along the
/// axis being worked on, [X|Y|Z]globmin, globmax, and the global positions of the boxes and planes.
/// A step that returns true has failed, why says why.
struct LatticeSetup
{
	const fdtd::grid::GlobalPlanes &planes;
	std::string &why;
	int here[3] = {}; // this rank's position in the rank grid
	// first global cell, PML included, of the rank at position idx along an axis
	static int rank_start(const std::vector<int> &cells, int idx) { int start = 0; for (int m = 0; m < idx; ++m) start += cells[m]; return start; }
	int i{};
	int j{};
	int k{};
	int m{};
	int s{};
	int xre{};
	int yre{};
	int zre{};
	int Xglobmin{};
	int Xglobmax{};
	int Yglobmin{};
	int Yglobmax{};
	int Zglobmin{};
	int Zglobmax{};
	int NRglobfarX{};
	int NRglobfarY{};
	int NRglobfarZ{};
	int NRlocalX{};
	int NRlocalY{};
	int NRlocalZ{};
	int NRglobneaX{};
	int NRglobneaY{};
	int NRglobneaZ{};
	int ARTglobneaY{};
	int ARTglobfarY{};

	int QSCAglobneaX{};
	int QSCAglobfarX{};
	int QSCAglobneaY{};
	int QSCAglobfarY{};
	int QSCAglobneaZ{};
	int QSCAglobfarZ{};

	int QTOTglobneaX{};
	int QTOTglobfarX{};
	int QTOTglobneaY{};
	int QTOTglobfarY{};
	int QTOTglobneaZ{};
	int QTOTglobfarZ{};

	int I_Animation{};
	int J_Animation{};
	int J_Animation2{};
	int J_Animation1{};
	int K_Animation{};

	int process{};
	int indices[3]{};
	double xnr{};
	double ynr{};
	double znr{};
	double a{};
	double b{};
	bool bTest1{};
	bool bTest2{};
	bool bExclTestX = false;
	bool bExclTestY = false;
	bool bExclTestZ = false;
	bool bNotAtAllCPMLX = false;
	bool bNotAtAllCPMLY = false;
	bool bNotAtAllCPMLZ = false;

	/// the domain decomposition: this rank's extent along each axis, its CPML ranges and type
	bool rank_extent();
	/// the TF/SF box, the ART planes, the flux boxes and the plot planes in global cell indices, given
	bool global_planes();
	/// the near-to-far box on this rank, its inner block and output offsets
	bool near_to_far();
	/// the TF/SF box with periodic boundaries; checks it fits the domain
	bool tfsf_box();
	/// the plot planes, ART planes, flux faces and TF/SF faces that lie on this rank
	bool local_planes();
	/// the plane-wave injection ranges of this rank
	bool plane_wave_faces();
	/// the scattering flux box, QSCA, on this rank
	bool scattering_box();
	/// the total flux box, QTOT, on this rank
	bool total_box();
};

bool LatticeSetup::rank_extent()
{
	subdomain = Subdomain();
	subdomain.rank = fdtd::grid::world_rank;
	fdtd::grid::axis_indices(fdtd::grid::world_rank, here);
	for (int a = 0; a < 3; ++a)
		subdomain.position[a] = here[a];

    //@gradient-decomposition
	//
	// DECOMPOSITION OF THE DOMAIN
	//
	// X axis
	//
	fdtd::grid::cells_per_rank_x.assign(num_mpi_processes_x, 0);
	if (equal_decomp_x)
	{ 
		xnr = lat.cells[0] / num_mpi_processes_x;
		xre = lat.cells[0] % num_mpi_processes_x;
		a = modf(xnr, &b); //in b the integer part of xnr, ynr, znr
		s = 0;

		for (m = 0; m < num_mpi_processes_x; ++m)
			fdtd::grid::cells_per_rank_x[m] = int(b);

		//redistribution of the remaining cells to the processes 
		while (xre > 0)
		{
			fdtd::grid::cells_per_rank_x[s] += 1;
			--xre;
			++s;
			if (s >= num_mpi_processes_x)
				s = 0;
		}
	}
	else
	{
		for (m = 0; m < num_mpi_processes_x; ++m)
		{
			fdtd::grid::cells_per_rank_x[m] = (m == num_mpi_processes_x - 1) ? (lat.cells[0] - arr_decomp_x[m]) 
												: (arr_decomp_x[m + 1] - arr_decomp_x[m]);
		}
	}

	//
	// Y axis
	//
	fdtd::grid::cells_per_rank_y.assign(num_mpi_processes_y, 0);
	if (equal_decomp_y)
	{ 
		ynr = lat.cells[1] / num_mpi_processes_y;
		yre = lat.cells[1] % num_mpi_processes_y;
		a = modf(ynr, &b);
		s = 0;

		for (m = 0; m < num_mpi_processes_y; ++m)
			fdtd::grid::cells_per_rank_y[m] = int(b);
	
		while (yre > 0)
		{
			fdtd::grid::cells_per_rank_y[s] += 1;
			--yre;
			++s;
			if (s >= num_mpi_processes_y)
				s = 0;
		}
	}
	else
	{
		for (m = 0; m < num_mpi_processes_y; ++m)
		{
			fdtd::grid::cells_per_rank_y[m] = (m == num_mpi_processes_y - 1) ? (lat.cells[1] - arr_decomp_y[m]) 
												: (arr_decomp_y[m + 1] - arr_decomp_y[m]);
		}
	}

	//
	// Z axis
	//
	fdtd::grid::cells_per_rank_z.assign(num_mpi_processes_z, 0);
	if (equal_decomp_z)
	{
		znr = lat.cells[2] / num_mpi_processes_z;
		zre = lat.cells[2] % num_mpi_processes_z;
		a = modf(znr, &b);
		s = 0;

		for (m = 0; m < num_mpi_processes_z; ++m)
			fdtd::grid::cells_per_rank_z[m] = int(b);

		while (zre > 0)
		{
			fdtd::grid::cells_per_rank_z[s] += 1;
			--zre;
			++s;
			if (s >= num_mpi_processes_z)
				s = 0;
		}
	}
	else
	{
		for (m = 0; m < num_mpi_processes_z; ++m)
		{
			fdtd::grid::cells_per_rank_z[m] = (m == num_mpi_processes_z - 1) ? (lat.cells[2] - arr_decomp_z[m]) 
												: (arr_decomp_z[m + 1] - arr_decomp_z[m]);
		}
	}
	//
	//
	// Start(Nea) and End(Far) of the physical domain in the global coordinate system, including CPML
	NRglobneaX = lat.cpml[0];
	NRglobneaY = lat.cpml[1];
	NRglobneaZ = lat.cpml[2];
	NRglobfarX = lat.cells[0] - lat.cpml[0];
	NRglobfarY = lat.cells[1] - lat.cpml[1];
	NRglobfarZ = lat.cells[2] - lat.cpml[2];

//
//
//
// Identifying the rank of the processes that DO NOT PARTICIPATE on the file output
// These processes have to be excluded from COMM_WORLD in order to create a communicator COMM_OUTPUT, etc. etc. etc.
// ...
// ...
// TYPES OF PROCESSES/SUBDOMAINS:
//   ONLY_CPML, They are ONLY_NEAR or ONLY_FAR,
//   NON_CPML, They are NORMAL
//   MIXED, They are mixed, such as NEAR/NORMAL, NORMAL/FAR or NEAR/NORMAL/FAR
//
//
//
	// this rank's extent along each axis, [globmin, globmax), and its CPML ranges
	i = here[0];
	j = here[1];
	k = here[2];
	process = fdtd::grid::world_rank;
	Xglobmin = rank_start(fdtd::grid::cells_per_rank_x, here[0]);
	NRlocalX = fdtd::grid::cells_per_rank_x[i];
	Xglobmax = Xglobmin + NRlocalX;
	Yglobmin = rank_start(fdtd::grid::cells_per_rank_y, here[1]);
	NRlocalY = fdtd::grid::cells_per_rank_y[j];
	Yglobmax = Yglobmin + NRlocalY;
	Zglobmin = rank_start(fdtd::grid::cells_per_rank_z, here[2]);
	NRlocalZ = fdtd::grid::cells_per_rank_z[k];
	Zglobmax = Zglobmin + NRlocalZ;
	const AxisCpml cx = axis_cpml_ranges(Xglobmin, NRlocalX, NRglobneaX, NRglobfarX, lat.cells[0]);
	const AxisCpml cy = axis_cpml_ranges(Yglobmin, NRlocalY, NRglobneaY, NRglobfarY, lat.cells[1]);
	const AxisCpml cz = axis_cpml_ranges(Zglobmin, NRlocalZ, NRglobneaZ, NRglobfarZ, lat.cells[2]);
	cx.apply(subdomain.i_cpml_min, subdomain.i_cpml_max, subdomain.i_nea_cpml_end, subdomain.i_far_cpml_start, subdomain.cpml_pad_nea_x,
			 subdomain.cpml_pad_far_x, subdomain.inner_dim_x, subdomain.inner_org_x, subdomain.which_type_of_process_x, subdomain.shift_output_x);
	cy.apply(subdomain.j_cpml_min, subdomain.j_cpml_max, subdomain.j_nea_cpml_end, subdomain.j_far_cpml_start, subdomain.cpml_pad_nea_y,
			 subdomain.cpml_pad_far_y, subdomain.inner_dim_y, subdomain.inner_org_y, subdomain.which_type_of_process_y, subdomain.shift_output_y);
	cz.apply(subdomain.k_cpml_min, subdomain.k_cpml_max, subdomain.k_nea_cpml_end, subdomain.k_far_cpml_start, subdomain.cpml_pad_nea_z,
			 subdomain.cpml_pad_far_z, subdomain.inner_dim_z, subdomain.inner_org_z, subdomain.which_type_of_process_z, subdomain.shift_output_z);
	bExclTestX = cx.only_cpml, bExclTestY = cy.only_cpml, bExclTestZ = cz.only_cpml;
	bNotAtAllCPMLX = cx.no_cpml, bNotAtAllCPMLY = cy.no_cpml, bNotAtAllCPMLZ = cz.no_cpml;

	// TO IDENTIFY THE RANKS OF ONLY CPML PROCESSES, WHICH WILL BE EXCLUDED FROM fdtd::grid::output_comm

	process = fdtd::grid::world_rank;
	{

		subdomain.not_cpml_at_all = false;
		subdomain.cpml_only = false;

		if (bExclTestX || bExclTestY || bExclTestZ)
		{

			subdomain.cpml_only = true;

		}

		if (bNotAtAllCPMLX && bNotAtAllCPMLY && bNotAtAllCPMLZ)

			subdomain.not_cpml_at_all = true;
	}

	return false;
}

bool LatticeSetup::global_planes()
{
	// the NRglobnea / NRglobfar variables, the CPML boundary above, are re-used for the TF/SF box
	NRglobneaX = planes.tfsf_nea[0];
	NRglobneaY = planes.tfsf_nea[1];
	NRglobneaZ = planes.tfsf_nea[2];
	NRglobfarX = planes.tfsf_far[0];
	NRglobfarY = planes.tfsf_far[1];
	NRglobfarZ = planes.tfsf_far[2];
	ARTglobneaY = planes.art_nea_y;
	ARTglobfarY = planes.art_far_y;
	QSCAglobneaX = planes.qsca_nea[0];
	QSCAglobfarX = planes.qsca_far[0];
	QTOTglobneaX = planes.qtot_nea[0];
	QTOTglobfarX = planes.qtot_far[0];
	QSCAglobneaY = planes.qsca_nea[1];
	QSCAglobfarY = planes.qsca_far[1];
	QTOTglobneaY = planes.qtot_nea[1];
	QTOTglobfarY = planes.qtot_far[1];
	QSCAglobneaZ = planes.qsca_nea[2];
	QSCAglobfarZ = planes.qsca_far[2];
	QTOTglobneaZ = planes.qtot_nea[2];
	QTOTglobfarZ = planes.qtot_far[2];
	I_Animation = planes.plot_x;
	J_Animation = planes.plot_y;
	K_Animation = planes.plot_z;
	J_Animation1 = planes.plot_y1;
	J_Animation2 = planes.plot_y2;

	if (fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("grid", "plot planes: x = %d, y = %d, z = %d, y1 = %d, y2 = %d\n", I_Animation, J_Animation, K_Animation, J_Animation1, J_Animation2);
	if (fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("grid", "flux boxes along y: scattering %d to %d, total %d to %d\n", QSCAglobneaY, QSCAglobfarY, QTOTglobneaY, QTOTglobfarY);

	return false;
}

bool LatticeSetup::near_to_far()
{

// for NEAR-to-Far plane extraction
	// this rank along x: [Xglobmin, Xglobmax)
	i = here[0];
	j = here[1];
	k = here[2];
	Xglobmin = rank_start(fdtd::grid::cells_per_rank_x, here[0]);

		NRlocalX = fdtd::grid::cells_per_rank_x[i];

		Xglobmax = Xglobmin + NRlocalX;

	{

				process = fdtd::grid::world_rank;

				if (Xglobmin <= QSCAglobneaX && Xglobmax > QSCAglobfarX)
				{
					// MIXED NEAR/NORMAL/FAR
					subdomain.inner_org_x_nf = QSCAglobneaX - Xglobmin;
				}
				else if (Xglobmin <= QSCAglobneaX && Xglobmax > QSCAglobneaX)
				{
					// MIXED NEAR/NORMAL
					subdomain.inner_org_x_nf = QSCAglobneaX - Xglobmin;
				}
				else if (Xglobmin > QSCAglobneaX && Xglobmax <= QSCAglobfarX)
				{
					// NON_CPML
					subdomain.inner_org_x_nf = 0;
				}
				else if (Xglobmin <= QSCAglobfarX && Xglobmax > QSCAglobfarX)
				{

					subdomain.inner_org_x_nf = 0;
				}

				// Referred to the global origin of the array of data
				subdomain.shift_output_x_nf = Xglobmin + subdomain.inner_org_x_nf - QSCAglobneaX;
				}

	// this rank along y: [Yglobmin, Yglobmax)
	i = here[0];
	j = here[1];
	k = here[2];
	Yglobmin = rank_start(fdtd::grid::cells_per_rank_y, here[1]);

		NRlocalY = fdtd::grid::cells_per_rank_y[j];

		Yglobmax = Yglobmin + NRlocalY;

	{

				process = fdtd::grid::world_rank;

				if (Yglobmin <= QSCAglobneaY && Yglobmax > QSCAglobfarY)
				{
					// MIXED NEAR/NORMAL/FAR
					subdomain.inner_org_y_nf = QSCAglobneaY - Yglobmin;
				}
				else if (Yglobmin <= QSCAglobneaY && Yglobmax > QSCAglobneaY)
				{
					// MIXED NEAR/NORMAL
					subdomain.inner_org_y_nf = QSCAglobneaY - Yglobmin;
				}
				else if (Yglobmin > QSCAglobneaY && Yglobmax <= QSCAglobfarY)
				{
					// NON_CPML
					subdomain.inner_org_y_nf = 0;
				}
				else if (Yglobmin <= QSCAglobfarY && Yglobmax > QSCAglobfarY)
				{
					// MIXED NORMAL/FAR
					subdomain.inner_org_y_nf = 0;
				}
				// Referred to the global origin of the array of data
				subdomain.shift_output_y_nf = Yglobmin + subdomain.inner_org_y_nf - QSCAglobneaY;
				}

	// this rank along z: [Zglobmin, Zglobmax)
	i = here[0];
	j = here[1];
	k = here[2];
	Zglobmin = rank_start(fdtd::grid::cells_per_rank_z, here[2]);

		NRlocalZ = fdtd::grid::cells_per_rank_z[k];

		Zglobmax = Zglobmin + NRlocalZ;

	{

				process = fdtd::grid::world_rank;

				if (Zglobmin <= QSCAglobneaZ && Zglobmax > QSCAglobfarZ)
				{
					// MIXED NEAR/NORMALE/FAR
					subdomain.inner_org_z_nf = QSCAglobneaZ - Zglobmin;
				}
				else if (Zglobmin <= QSCAglobneaZ && Zglobmax > QSCAglobneaZ)
				{
					// MIXED NEAR/NORMALE
					subdomain.inner_org_z_nf = QSCAglobneaZ - Zglobmin;
				}
				else if (Zglobmin > QSCAglobneaZ && Zglobmax <= QSCAglobfarZ)
				{
					// NON_CPML
					subdomain.inner_org_z_nf = 0;
				}
				else if (Zglobmin <= QSCAglobfarZ && Zglobmax > QSCAglobfarZ)
				{
					// MIXED NORMALE/FAR
					subdomain.inner_org_z_nf = 0;
				}
				// Referred to the global origin of the array of data
				subdomain.shift_output_z_nf = Zglobmin + subdomain.inner_org_z_nf - QSCAglobneaZ;
				}

	return false;
}

bool LatticeSetup::tfsf_box()
{
	if (lat.periodic_x == 0)
	{
		Xglobmin = 0;
		for (i = 0; i < num_mpi_processes_x; ++i)
		{

			NRlocalX = fdtd::grid::cells_per_rank_x[i];

			Xglobmax = Xglobmin + NRlocalX;

			bTest1 = (NRglobneaX == Xglobmin) || (NRglobneaX == Xglobmax);
			bTest2 = (NRglobfarX == Xglobmin) || (NRglobfarX == Xglobmax);

			if (bTest1 || bTest2)
			{
				// rank_starts() (sim_config) keeps rank boundaries off the TF/SF planes; reaching this is a bug
				if (fdtd::grid::world_rank == ROOTPROCESS)
					fdtd::io::log("grid", "internal error: a %s rank boundary lies on the TF/SF plane\n", "x");
				why = "a rank boundary lies on the TF/SF plane";
				return true;
			}

			Xglobmin = Xglobmax;
		}
	}
	Yglobmin = 0;
	for (j = 0; j < num_mpi_processes_y; ++j)
	{

		NRlocalY = fdtd::grid::cells_per_rank_y[j];

		Yglobmax = Yglobmin + NRlocalY;

		bTest1 = (NRglobneaY == Yglobmin) || (NRglobneaY == Yglobmax);
		bTest2 = (NRglobfarY == Yglobmin) || (NRglobfarY == Yglobmax);

		if (bTest1 || bTest2)
		{
			// rank_starts() (sim_config) keeps rank boundaries off the TF/SF planes; reaching this is a bug
			if (fdtd::grid::world_rank == ROOTPROCESS)
				fdtd::io::log("grid", "internal error: a %s rank boundary lies on the TF/SF plane\n", "y");
			why = "a rank boundary lies on the TF/SF plane";
			return true;
		}

		Yglobmin = Yglobmax;
	}

	if (lat.periodic_z == 0)
	{
		Zglobmin = 0;
		for (k = 0; k < num_mpi_processes_z; ++k)
		{

			NRlocalZ = fdtd::grid::cells_per_rank_z[k];

			Zglobmax = Zglobmin + NRlocalZ;

			bTest1 = (NRglobneaZ == Zglobmin) || (NRglobneaZ == Zglobmax);
			bTest2 = (NRglobfarZ == Zglobmin) || (NRglobfarZ == Zglobmax);

			if (bTest1 || bTest2)
			{
				// rank_starts() (sim_config) keeps rank boundaries off the TF/SF planes; reaching this is a bug
				if (fdtd::grid::world_rank == ROOTPROCESS)
					fdtd::io::log("grid", "internal error: a %s rank boundary lies on the TF/SF plane\n", "z");
				why = "a rank boundary lies on the TF/SF plane";
				return true;
			}

			Zglobmin = Zglobmax;
		}
	}
	return false;
}

bool LatticeSetup::local_planes()
{
	// this rank along x: [Xglobmin, Xglobmax)
	i = here[0];
	j = here[1];
	k = here[2];
	Xglobmin = rank_start(fdtd::grid::cells_per_rank_x, here[0]);

		NRlocalX = fdtd::grid::cells_per_rank_x[i];

		Xglobmax = Xglobmin + NRlocalX;

	{

				process = fdtd::grid::world_rank;

				if ((subdomain.pw_nea_x =
					((Xglobmin <= NRglobneaX) && (NRglobneaX < Xglobmax))))
				{
					subdomain.ipw_1 = NRglobneaX - Xglobmin;
				}

				if ((subdomain.pw_far_x =
					((Xglobmin <= NRglobfarX) && (NRglobfarX < Xglobmax))))
				{
					subdomain.ipw_2 = NRglobfarX - Xglobmin;
				}

				if ((subdomain.qsca_nea_x =
					((Xglobmin <= QSCAglobneaX) && (QSCAglobneaX < Xglobmax))))
				{
					subdomain.i_qsca_nea = int(QSCAglobneaX - Xglobmin - subdomain.inner_org_x);
				}

				if ((subdomain.qsca_far_x =
					((Xglobmin <= QSCAglobfarX) && (QSCAglobfarX < Xglobmax))))
				{
					subdomain.i_qsca_far = int(QSCAglobfarX - Xglobmin - subdomain.inner_org_x);
				}

				if ((subdomain.qtot_nea_x =
					((Xglobmin <= QTOTglobneaX) && (QTOTglobneaX < Xglobmax))))
				{
					subdomain.i_qtot_nea = int(QTOTglobneaX - Xglobmin - subdomain.inner_org_x);
				}

				if ((subdomain.qtot_far_x =
					((Xglobmin <= QTOTglobfarX) && (QTOTglobfarX < Xglobmax))))
				{
					subdomain.i_qtot_far = int(QTOTglobfarX - Xglobmin - subdomain.inner_org_x);
				}

				if ((subdomain.animation_yz =
					((Xglobmin <= I_Animation) && (I_Animation < Xglobmax))))
				{
					subdomain.i_animation = int(I_Animation - Xglobmin - subdomain.inner_org_x);
				}
				}

	// this rank along y: [Yglobmin, Yglobmax)
	i = here[0];
	j = here[1];
	k = here[2];
	Yglobmin = rank_start(fdtd::grid::cells_per_rank_y, here[1]);

		NRlocalY = fdtd::grid::cells_per_rank_y[j];

		Yglobmax = Yglobmin + NRlocalY;

	{

				process = fdtd::grid::world_rank;

				if ((subdomain.pw_nea_y =
					((Yglobmin <= NRglobneaY) && (NRglobneaY < Yglobmax))))
				{
					subdomain.jpw_1 = NRglobneaY - Yglobmin;
				}

				if ((subdomain.pw_far_y =
					((Yglobmin <= NRglobfarY) && (NRglobfarY < Yglobmax))))
				{
					subdomain.jpw_2 = NRglobfarY - Yglobmin;
				}

				if ((subdomain.qsca_nea_y =
					((Yglobmin <= QSCAglobneaY) && (QSCAglobneaY < Yglobmax))))
				{
					subdomain.j_qsca_nea = int(QSCAglobneaY - Yglobmin - subdomain.inner_org_y);
				}

				if ((subdomain.qsca_far_y =
					((Yglobmin <= QSCAglobfarY) && (QSCAglobfarY < Yglobmax))))
				{
					subdomain.j_qsca_far = int(QSCAglobfarY - Yglobmin - subdomain.inner_org_y);
				}

				if ((subdomain.qtot_nea_y =
					((Yglobmin <= QTOTglobneaY) && (QTOTglobneaY < Yglobmax))))
				{
					subdomain.j_qtot_nea = int(QTOTglobneaY - Yglobmin - subdomain.inner_org_y);
				}

				if ((subdomain.qtot_far_y =
					((Yglobmin <= QTOTglobfarY) && (QTOTglobfarY < Yglobmax))))
				{
					subdomain.j_qtot_far = int(QTOTglobfarY - Yglobmin - subdomain.inner_org_y);
				}

				if ((subdomain.art_nea_y =
					((Yglobmin <= ARTglobneaY) && (ARTglobneaY < Yglobmax))))
				{
					subdomain.j_art_1 = ARTglobneaY - Yglobmin - subdomain.inner_org_y;
				}

				if ((subdomain.art_far_y =
					((Yglobmin <= ARTglobfarY) && (ARTglobfarY < Yglobmax))))
				{
					subdomain.j_art_2 = ARTglobfarY - Yglobmin - subdomain.inner_org_y;
				}

				if ((subdomain.animation_xz =
					((Yglobmin <= J_Animation) && (J_Animation < Yglobmax))))
				{
					subdomain.j_animation = int(J_Animation - Yglobmin - subdomain.inner_org_y);
					if (planes.plot_y_not_on_rank_start && subdomain.j_animation == 0)
					{
						why = "graphene_optim: the animation plane falls on the first cell of a rank";
						return true;
					}
				}

				if ((subdomain.animation_xz_2 =
					((Yglobmin <= J_Animation2) && (J_Animation2 < Yglobmax))))
				{
					subdomain.j_animation_2 = int(J_Animation2 - Yglobmin - subdomain.inner_org_y);
				}

				if ((subdomain.animation_xz_1 =
					((Yglobmin <= J_Animation1) && (J_Animation1 < Yglobmax))))
				{
					subdomain.j_animation_1 = int(J_Animation1 - Yglobmin - subdomain.inner_org_y);
				}
				}

	// this rank along z: [Zglobmin, Zglobmax)
	i = here[0];
	j = here[1];
	k = here[2];
	Zglobmin = rank_start(fdtd::grid::cells_per_rank_z, here[2]);

		NRlocalZ = fdtd::grid::cells_per_rank_z[k];

		Zglobmax = Zglobmin + NRlocalZ;

	{

				process = fdtd::grid::world_rank;

				if ((subdomain.pw_nea_z =
					((Zglobmin <= NRglobneaZ) && (NRglobneaZ < Zglobmax))))
				{
					subdomain.kpw_1 = NRglobneaZ - Zglobmin;
				}

				if ((subdomain.pw_far_z =
					((Zglobmin <= NRglobfarZ) && (NRglobfarZ < Zglobmax))))
				{
					subdomain.kpw_2 = NRglobfarZ - Zglobmin;
				}

				if ((subdomain.qsca_nea_z =
					((Zglobmin <= QSCAglobneaZ) && (QSCAglobneaZ < Zglobmax))))
				{
					subdomain.k_qsca_nea = int(QSCAglobneaZ - Zglobmin - subdomain.inner_org_z);
				}

				if ((subdomain.qsca_far_z =
					((Zglobmin <= QSCAglobfarZ) && (QSCAglobfarZ < Zglobmax))))
				{
					subdomain.k_qsca_far = int(QSCAglobfarZ - Zglobmin - subdomain.inner_org_z);
				}

				if ((subdomain.qtot_nea_z =
					((Zglobmin <= QTOTglobneaZ) && (QTOTglobneaZ < Zglobmax))))
				{
					subdomain.k_qtot_nea = int(QTOTglobneaZ - Zglobmin - subdomain.inner_org_z);
				}

				if ((subdomain.qtot_far_z =
					((Zglobmin <= QTOTglobfarZ) && (QTOTglobfarZ < Zglobmax))))
				{
					subdomain.k_qtot_far = int(QTOTglobfarZ - Zglobmin - subdomain.inner_org_z);
				}

				if ((subdomain.animation_xy =
					((Zglobmin <= K_Animation) && (K_Animation < Zglobmax))))
				{
					subdomain.k_animation = int(K_Animation - Zglobmin - subdomain.inner_org_z);
				}
				}

	return false;
}

bool LatticeSetup::plane_wave_faces()
{
	process = fdtd::grid::world_rank;
	{

		fdtd::grid::axis_indices(process, indices);

		if (subdomain.pw_nea_x || subdomain.pw_far_x)
		{

			Yglobmin = 0;
			for (m = 0; m < indices[1]; ++m)
				Yglobmin += fdtd::grid::cells_per_rank_y[m];
			NRlocalY = fdtd::grid::cells_per_rank_y[m];
			Yglobmax = Yglobmin + NRlocalY;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.jpw_min = std::max(std::max(NRglobneaY, Yglobmin) - Yglobmin, 0);
			subdomain.jpw_max = std::min(NRglobfarY, Yglobmax) - Yglobmin;
			subdomain.jpw_max = std::max(subdomain.jpw_max, subdomain.jpw_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(NRglobneaY, Yglobmin) > std::min(NRglobfarY, Yglobmax))
				subdomain.pw_nea_x = subdomain.pw_far_x = false;

			Zglobmin = 0;
			for (m = 0; m < indices[2]; ++m)
				Zglobmin += fdtd::grid::cells_per_rank_z[m];
			NRlocalZ = fdtd::grid::cells_per_rank_z[m];
			Zglobmax = Zglobmin + NRlocalZ;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.kpw_min = std::max(std::max(NRglobneaZ, Zglobmin) - Zglobmin, 0);
			subdomain.kpw_max = std::min(NRglobfarZ, Zglobmax) - Zglobmin;
			subdomain.kpw_max = std::max(subdomain.kpw_max, subdomain.kpw_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(NRglobneaZ, Zglobmin) > std::min(NRglobfarZ, Zglobmax))
				subdomain.pw_nea_x = subdomain.pw_far_x = false;
		}

		if (subdomain.pw_nea_y || subdomain.pw_far_y)
		{

			Xglobmin = 0;
			for (m = 0; m < indices[0]; ++m)
				Xglobmin += fdtd::grid::cells_per_rank_x[m];
			NRlocalX = fdtd::grid::cells_per_rank_x[m];
			Xglobmax = Xglobmin + NRlocalX;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.ipw_min = std::max(std::max(NRglobneaX, Xglobmin) - Xglobmin, 0);
			subdomain.ipw_max = std::min(NRglobfarX, Xglobmax) - Xglobmin;
			subdomain.ipw_max = std::max(subdomain.ipw_max, subdomain.ipw_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(NRglobneaX, Xglobmin) > std::min(NRglobfarX, Xglobmax))
				subdomain.pw_nea_y = subdomain.pw_far_y = false;

			Zglobmin = 0;
			for (m = 0; m < indices[2]; ++m)
				Zglobmin += fdtd::grid::cells_per_rank_z[m];
			NRlocalZ = fdtd::grid::cells_per_rank_z[m];
			Zglobmax = Zglobmin + NRlocalZ;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.kpw_min = std::max(std::max(NRglobneaZ, Zglobmin) - Zglobmin, 0);
			subdomain.kpw_max = std::min(NRglobfarZ, Zglobmax) - Zglobmin;
			subdomain.kpw_max = std::max(subdomain.kpw_max, subdomain.kpw_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(NRglobneaZ, Zglobmin) > std::min(NRglobfarZ, Zglobmax))
				subdomain.pw_nea_y = subdomain.pw_far_y = false;
		}

		if (subdomain.pw_nea_z || subdomain.pw_far_z)
		{

			Xglobmin = 0;
			for (m = 0; m < indices[0]; ++m)
				Xglobmin += fdtd::grid::cells_per_rank_x[m];
			NRlocalX = fdtd::grid::cells_per_rank_x[m];
			Xglobmax = Xglobmin + NRlocalX;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.ipw_min = std::max(std::max(NRglobneaX, Xglobmin) - Xglobmin, 0);
			subdomain.ipw_max = std::min(NRglobfarX, Xglobmax) - Xglobmin;
			subdomain.ipw_max = std::max(subdomain.ipw_max, subdomain.ipw_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(NRglobneaX, Xglobmin) > std::min(NRglobfarX, Xglobmax))
				subdomain.pw_nea_z = subdomain.pw_far_z = false;

			Yglobmin = 0;
			for (m = 0; m < indices[1]; ++m)
				Yglobmin += fdtd::grid::cells_per_rank_y[m];
			NRlocalY = fdtd::grid::cells_per_rank_y[m];
			Yglobmax = Yglobmin + NRlocalY;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.jpw_min = std::max(std::max(NRglobneaY, Yglobmin) - Yglobmin, 0);
			subdomain.jpw_max = std::min(NRglobfarY, Yglobmax) - Yglobmin;
			subdomain.jpw_max = std::max(subdomain.jpw_max, subdomain.jpw_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(NRglobneaY, Yglobmin) > std::min(NRglobfarY, Yglobmax))
				subdomain.pw_nea_z = subdomain.pw_far_z = false;
		}

		Yglobmin = 0;
		for (m = 0; m < indices[1]; ++m)
			Yglobmin += fdtd::grid::cells_per_rank_y[m];
		NRlocalY = fdtd::grid::cells_per_rank_y[m];
		Yglobmax = Yglobmin + NRlocalY;
		subdomain.jpw_pad = Yglobmin;
	}
	return false;
}

bool LatticeSetup::scattering_box()
{
	// QSCA
	process = fdtd::grid::world_rank;
	{

		fdtd::grid::axis_indices(process, indices);

		if (subdomain.qsca_nea_x || subdomain.qsca_far_x)
		{

			Yglobmin = 0;
			for (m = 0; m < indices[1]; ++m)
				Yglobmin += fdtd::grid::cells_per_rank_y[m];
			NRlocalY = fdtd::grid::cells_per_rank_y[m];
			Yglobmax = Yglobmin + NRlocalY;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.j_qsca_min = std::max(std::max(QSCAglobneaY, Yglobmin) - Yglobmin - subdomain.inner_org_y, 0);
			subdomain.j_qsca_max = std::min(QSCAglobfarY, Yglobmax) - Yglobmin - subdomain.inner_org_y;
			subdomain.j_qsca_max = std::max(subdomain.j_qsca_max, subdomain.j_qsca_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QSCAglobneaY, Yglobmin) > std::min(QSCAglobfarY, Yglobmax))
				subdomain.qsca_nea_x = subdomain.qsca_far_x = false;

			Zglobmin = 0;
			for (m = 0; m < indices[2]; ++m)
				Zglobmin += fdtd::grid::cells_per_rank_z[m];
			NRlocalZ = fdtd::grid::cells_per_rank_z[m];
			Zglobmax = Zglobmin + NRlocalZ;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.k_qsca_min = std::max(std::max(QSCAglobneaZ, Zglobmin) - Zglobmin - subdomain.inner_org_z, 0);
			subdomain.k_qsca_max = std::min(QSCAglobfarZ, Zglobmax) - Zglobmin - subdomain.inner_org_z;
			subdomain.k_qsca_max = std::max(subdomain.k_qsca_max, subdomain.k_qsca_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QSCAglobneaZ, Zglobmin) > std::min(QSCAglobfarZ, Zglobmax))
				subdomain.qsca_nea_x = subdomain.qsca_far_x = false;
		}

		if (subdomain.qsca_nea_y || subdomain.qsca_far_y)
		{

			Xglobmin = 0;
			for (m = 0; m < indices[0]; ++m)
				Xglobmin += fdtd::grid::cells_per_rank_x[m];
			NRlocalX = fdtd::grid::cells_per_rank_x[m];
			Xglobmax = Xglobmin + NRlocalX;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.i_qsca_min = std::max(std::max(QSCAglobneaX, Xglobmin) - Xglobmin - subdomain.inner_org_x, 0);
			subdomain.i_qsca_max = std::min(QSCAglobfarX, Xglobmax) - Xglobmin - subdomain.inner_org_x;
			subdomain.i_qsca_max = std::max(subdomain.i_qsca_max, subdomain.i_qsca_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QSCAglobneaX, Xglobmin) > std::min(QSCAglobfarX, Xglobmax))
				subdomain.qsca_nea_y = subdomain.qsca_far_y = false;

			Zglobmin = 0;
			for (m = 0; m < indices[2]; ++m)
				Zglobmin += fdtd::grid::cells_per_rank_z[m];
			NRlocalZ = fdtd::grid::cells_per_rank_z[m];
			Zglobmax = Zglobmin + NRlocalZ;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.k_qsca_min = std::max(std::max(QSCAglobneaZ, Zglobmin) - Zglobmin - subdomain.inner_org_z, 0);
			subdomain.k_qsca_max = std::min(QSCAglobfarZ, Zglobmax) - Zglobmin - subdomain.inner_org_z;
			subdomain.k_qsca_max = std::max(subdomain.k_qsca_max, subdomain.k_qsca_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QSCAglobneaZ, Zglobmin) > std::min(QSCAglobfarZ, Zglobmax))
				subdomain.qsca_nea_y = subdomain.qsca_far_y = false;
		}

		if (subdomain.qsca_nea_z || subdomain.qsca_far_z)
		{

			Xglobmin = 0;
			for (m = 0; m < indices[0]; ++m)
				Xglobmin += fdtd::grid::cells_per_rank_x[m];
			NRlocalX = fdtd::grid::cells_per_rank_x[m];
			Xglobmax = Xglobmin + NRlocalX;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.i_qsca_min = std::max(std::max(QSCAglobneaX, Xglobmin) - Xglobmin - subdomain.inner_org_x, 0);
			subdomain.i_qsca_max = std::min(QSCAglobfarX, Xglobmax) - Xglobmin - subdomain.inner_org_x;
			subdomain.i_qsca_max = std::max(subdomain.i_qsca_max, subdomain.i_qsca_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QSCAglobneaX, Xglobmin) > std::min(QSCAglobfarX, Xglobmax))
				subdomain.qsca_nea_z = subdomain.qsca_far_z = false;

			Yglobmin = 0;
			for (m = 0; m < indices[1]; ++m)
				Yglobmin += fdtd::grid::cells_per_rank_y[m];
			NRlocalY = fdtd::grid::cells_per_rank_y[m];
			Yglobmax = Yglobmin + NRlocalY;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.j_qsca_min = std::max(std::max(QSCAglobneaY, Yglobmin) - Yglobmin - subdomain.inner_org_y, 0);
			subdomain.j_qsca_max = std::min(QSCAglobfarY, Yglobmax) - Yglobmin - subdomain.inner_org_y;
			subdomain.j_qsca_max = std::max(subdomain.j_qsca_max, subdomain.j_qsca_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QSCAglobneaY, Yglobmin) > std::min(QSCAglobfarY, Yglobmax))
				subdomain.qsca_nea_z = subdomain.qsca_far_z = false;
		}

		Yglobmin = 0;
		for (m = 0; m < indices[1]; ++m)
			Yglobmin += fdtd::grid::cells_per_rank_y[m];
		NRlocalY = fdtd::grid::cells_per_rank_y[m];
		Yglobmax = Yglobmin + NRlocalY;

		//if(subdomain.animation_yz && (  (subdomain.j_qsca_min != subdomain.j_qsca_max) &&
		//  subdomain.k_qsca_min != subdomain.k_qsca_max)) // && !subdomain.qsca_nea_y && !subdomain.qsca_nea_z && !subdomain.qsca_far_x && !subdomain.qsca_far_y && !subdomain.qsca_far_z)
	}

	return false;
}

bool LatticeSetup::total_box()
{
	//for NEAR-to-FAR

	// QTOT
	process = fdtd::grid::world_rank;
	{

		fdtd::grid::axis_indices(process, indices);

		if (subdomain.qtot_nea_x || subdomain.qtot_far_x)
		{

			Yglobmin = 0;
			for (m = 0; m < indices[1]; ++m)
				Yglobmin += fdtd::grid::cells_per_rank_y[m];
			NRlocalY = fdtd::grid::cells_per_rank_y[m];
			Yglobmax = Yglobmin + NRlocalY;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.j_qtot_min = std::max(std::max(QTOTglobneaY, Yglobmin) - Yglobmin - subdomain.inner_org_y, 0);
			subdomain.j_qtot_max = std::min(QTOTglobfarY, Yglobmax) - Yglobmin - subdomain.inner_org_y;
			subdomain.j_qtot_max = std::max(subdomain.j_qtot_max, subdomain.j_qtot_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QTOTglobneaY, Yglobmin) > std::min(QTOTglobfarY, Yglobmax))
				subdomain.qtot_nea_x = subdomain.qtot_far_x = false;

			Zglobmin = 0;
			for (m = 0; m < indices[2]; ++m)
				Zglobmin += fdtd::grid::cells_per_rank_z[m];
			NRlocalZ = fdtd::grid::cells_per_rank_z[m];
			Zglobmax = Zglobmin + NRlocalZ;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.k_qtot_min = std::max(std::max(QTOTglobneaZ, Zglobmin) - Zglobmin - subdomain.inner_org_z, 0);
			subdomain.k_qtot_max = std::min(QTOTglobfarZ, Zglobmax) - Zglobmin - subdomain.inner_org_z;
			subdomain.k_qtot_max = std::max(subdomain.k_qtot_max, subdomain.k_qtot_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QTOTglobneaZ, Zglobmin) > std::min(QTOTglobfarZ, Zglobmax))
				subdomain.qtot_nea_x = subdomain.qtot_far_x = false;
		}

		if (subdomain.qtot_nea_y || subdomain.qtot_far_y)
		{

			Xglobmin = 0;
			for (m = 0; m < indices[0]; ++m)
				Xglobmin += fdtd::grid::cells_per_rank_x[m];
			NRlocalX = fdtd::grid::cells_per_rank_x[m];
			Xglobmax = Xglobmin + NRlocalX;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.i_qtot_min = std::max(std::max(QTOTglobneaX, Xglobmin) - Xglobmin - subdomain.inner_org_x, 0);
			subdomain.i_qtot_max = std::min(QTOTglobfarX, Xglobmax) - Xglobmin - subdomain.inner_org_x;
			subdomain.i_qtot_max = std::max(subdomain.i_qtot_max, subdomain.i_qtot_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QTOTglobneaX, Xglobmin) > std::min(QTOTglobfarX, Xglobmax))
				subdomain.qtot_nea_y = subdomain.qtot_far_y = false;

			Zglobmin = 0;
			for (m = 0; m < indices[2]; ++m)
				Zglobmin += fdtd::grid::cells_per_rank_z[m];
			NRlocalZ = fdtd::grid::cells_per_rank_z[m];
			Zglobmax = Zglobmin + NRlocalZ;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.k_qtot_min = std::max(std::max(QTOTglobneaZ, Zglobmin) - Zglobmin - subdomain.inner_org_z, 0);
			subdomain.k_qtot_max = std::min(QTOTglobfarZ, Zglobmax) - Zglobmin - subdomain.inner_org_z;
			subdomain.k_qtot_max = std::max(subdomain.k_qtot_max, subdomain.k_qtot_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QTOTglobneaZ, Zglobmin) > std::min(QTOTglobfarZ, Zglobmax))
				subdomain.qtot_nea_y = subdomain.qtot_far_y = false;
		}

		if (subdomain.qtot_nea_z || subdomain.qtot_far_z)
		{

			Xglobmin = 0;
			for (m = 0; m < indices[0]; ++m)
				Xglobmin += fdtd::grid::cells_per_rank_x[m];
			NRlocalX = fdtd::grid::cells_per_rank_x[m];
			Xglobmax = Xglobmin + NRlocalX;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.i_qtot_min = std::max(std::max(QTOTglobneaX, Xglobmin) - Xglobmin - subdomain.inner_org_x, 0);
			subdomain.i_qtot_max = std::min(QTOTglobfarX, Xglobmax) - Xglobmin - subdomain.inner_org_x;
			subdomain.i_qtot_max = std::max(subdomain.i_qtot_max, subdomain.i_qtot_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QTOTglobneaX, Xglobmin) > std::min(QTOTglobfarX, Xglobmax))
				subdomain.qtot_nea_z = subdomain.qtot_far_z = false;

			Yglobmin = 0;
			for (m = 0; m < indices[1]; ++m)
				Yglobmin += fdtd::grid::cells_per_rank_y[m];
			NRlocalY = fdtd::grid::cells_per_rank_y[m];
			Yglobmax = Yglobmin + NRlocalY;

			// overlap of [near, far) with this rank, a rank may hold both planes, e.g. a single rank
			subdomain.j_qtot_min = std::max(std::max(QTOTglobneaY, Yglobmin) - Yglobmin - subdomain.inner_org_y, 0);
			subdomain.j_qtot_max = std::min(QTOTglobfarY, Yglobmax) - Yglobmin - subdomain.inner_org_y;
			subdomain.j_qtot_max = std::max(subdomain.j_qtot_max, subdomain.j_qtot_min);
			// the face only exists on this rank if the rank also overlaps it transversally
			if (std::max(QTOTglobneaY, Yglobmin) > std::min(QTOTglobfarY, Yglobmax))
				subdomain.qtot_nea_z = subdomain.qtot_far_z = false;
		}

		Yglobmin = 0;
		for (m = 0; m < indices[1]; ++m)
			Yglobmin += fdtd::grid::cells_per_rank_y[m];
		NRlocalY = fdtd::grid::cells_per_rank_y[m];
		Yglobmax = Yglobmin + NRlocalY;
	}

	return false;
}

} // namespace

bool fdtd::grid::set_lattice_data(const GlobalPlanes& planes, std::string& why)
{
	LatticeSetup s{planes, why};
	return s.rank_extent() || s.global_planes() || s.near_to_far() || s.tfsf_box() || s.local_planes() || s.plane_wave_faces() || s.scattering_box() || s.total_box();
}

/**
 * @brief given the rank, the indices i, j, k are found. These
 * are the spatial positioning of the subdomain
 * 
 * @param fdtd::grid::world_rank rank of the current process
 * @param subdomain.position [out] This is where the determined indices are stored
 * 
 * @details This mapping is used for the sequential numbering of the MPI
 * processes.
 * N = i + x_procs_size * j + x_procs_size * y_procs_size * k \n
 * i = 0, ... , x_procs_size - 1 \n
 * j = 0, ... , y_procs_size - 1 \n
 * k = 0, ... , z_procs_size - 1 \n
 */
void fdtd::grid::axis_indices(const int &rank, int pos[3])
{
	int NX = num_mpi_processes_x;
	int NY = num_mpi_processes_y;
	int NZ = num_mpi_processes_z;
	int m;
	int N = NX * NY;
	double int_part;
	double frac_part;
	double u;

	u = double(rank) / double(N);

	for (m = NZ; m > 0; --m)
	{
		if (fabs(double(m) - u) < 1.E-08)
		{
			u = floor(u + 0.5);
			break;
		}
	}

	frac_part = modf(u, &int_part);

	pos[2] = int(int_part);

	u = frac_part * double(NY);

	for (m = NY; m > 0; --m)
	{
		if (fabs(double(m) - u) < 1.E-08)
		{
			u = floor(u + 0.5);
			break;
		}
	}

	frac_part = modf(u, &int_part);

	pos[1] = int(int_part);

	u = frac_part * double(NX);

	for (m = NX; m > 0; --m)
	{
		if (fabs(double(m) - u) < 1.E-08)
		{
			u = floor(u + 0.5);
			break;
		}
	}

	frac_part = modf(u, &int_part);

	pos[0] = int(int_part);

	return;
}

/**
 * @brief Uses Axis indices to determine the rank of neighbouring cells in
 * each direction
 * 
 * @param subdomain.position Array with i, j, k of current subdomain 
 * @param subdomain.neighbours_x Array to hold rank of each neighbour in X direction
 * @param subdomain.neighbours_y Array to hold rank of each neighbour in Y direction
 * @param subdomain.neighbours_z Array to hold rank of each neighbour in Z direction
 */
void fdtd::grid::neighbour_ranks(int *pos, int *nx, int *ny, int *nz)
{
	int i_less = pos[0] - 1;
	int i_more = pos[0] + 1;
	int j_less = pos[1] - 1;
	int j_more = pos[1] + 1;
	int k_less = pos[2] - 1;
	int k_more = pos[2] + 1;

	if (i_less >= 0)
		nx[0] = fdtd::grid::process_rank(i_less, pos[1], pos[2]);
	else
	{
		if (lat.periodic_x)
			nx[0] = fdtd::grid::process_rank(num_mpi_processes_x - 1, pos[1], pos[2]);
		else
			nx[0] = MPI_PROC_NULL;
	}
	if (i_more < num_mpi_processes_x)
		nx[1] = fdtd::grid::process_rank(i_more, pos[1], pos[2]);
	else
	{
		if (lat.periodic_x)
			nx[1] = fdtd::grid::process_rank(0, pos[1], pos[2]);
		else
			nx[1] = MPI_PROC_NULL;
	}
	if (j_less >= 0)
		ny[0] = fdtd::grid::process_rank(pos[0], j_less, pos[2]);
	else
		ny[0] = MPI_PROC_NULL;

	if (j_more < num_mpi_processes_y)
		ny[1] = fdtd::grid::process_rank(pos[0], j_more, pos[2]);
	else
		ny[1] = MPI_PROC_NULL;

	if (k_less >= 0)
		nz[0] = fdtd::grid::process_rank(pos[0], pos[1], k_less);
	else
	{
		if (lat.periodic_z)
			nz[0] = fdtd::grid::process_rank(pos[0], pos[1], num_mpi_processes_z - 1);
		else
			nz[0] = MPI_PROC_NULL;
	}
	if (k_more < num_mpi_processes_z)
		nz[1] = fdtd::grid::process_rank(pos[0], pos[1], k_more);
	else
	{
		if (lat.periodic_z)
			nz[1] = fdtd::grid::process_rank(pos[0], pos[1], 0);
		else
			nz[1] = MPI_PROC_NULL;
	}
	return;
}

/**
 * @brief Get the Neighbours Ranks N L object
 * 
 * @param subdomain.position Array with i, j, k of current subdomain 
 * @param subdomain.neighbours_nl_x Array to hold rank of neighbours in X direction
 * @param subdomain.neighbours_nl_y Array to hold rank of neighbours in Y direction
 * @param subdomain.neighbours_nl_z Array to hold rank of neighbours in Z direction
 * 
 * @note This was created to implement the sendrecv only in a box contained in
 * the simulation domain where the box encloses the nanostructure - this was
 * useful for nl simulations
 */
void fdtd::grid::neighbour_ranks_nl(int *pos, int *nx, int *ny, int *nz)
{
	// the plane-wave faces of every rank, each rank knows only its own
	const int mine[6] = {subdomain.pw_nea_x, subdomain.pw_far_x, subdomain.pw_nea_y, subdomain.pw_far_y, subdomain.pw_nea_z, subdomain.pw_far_z};
	std::vector<int> faces(6 * fdtd::grid::world_size);
	MPI_Allgather(mine, 6, MPI_INT, faces.data(), 6, MPI_INT, MPI_COMM_WORLD);
	auto face = [&](int which, int proc) { return faces[6 * proc + which] != 0; };
	int i;
	int j;
	int k;
	int process;
	// the ranks holding the TF/SF faces, no face found: no limit on that side
	int ind_PWNEAX = -1;
	int ind_PWNEAY = -1;
	int ind_PWNEAZ = -1;
	int ind_PWFARX = num_mpi_processes_x;
	int ind_PWFARY = num_mpi_processes_y;
	int ind_PWFARZ = num_mpi_processes_z;
	int i_less = pos[0] - 1;
	int i_more = pos[0] + 1;
	int j_less = pos[1] - 1;
	int j_more = pos[1] + 1;
	int k_less = pos[2] - 1;
	int k_more = pos[2] + 1;

	if (!lat.periodic_x && !lat.periodic_z)
	{
		j = int(0.5 * num_mpi_processes_y);
		k = int(0.5 * num_mpi_processes_z);

		for (i = 0; i < num_mpi_processes_x; ++i)
		{
			process = fdtd::grid::process_rank(i, j, k);
			if (face(0, process))
				ind_PWNEAX = i;
			if (face(1, process))
				ind_PWFARX = i;
		}

		i = int(0.5 * num_mpi_processes_x);
		k = int(0.5 * num_mpi_processes_z);

		for (j = 0; j < num_mpi_processes_y; ++j)
		{
			process = fdtd::grid::process_rank(i, j, k);
			if (face(2, process))
				ind_PWNEAY = j;
			if (face(3, process))
				ind_PWFARY = j;
		}

		i = int(0.5 * num_mpi_processes_x);
		j = int(0.5 * num_mpi_processes_y);

		for (k = 0; k < num_mpi_processes_z; ++k)
		{
			process = fdtd::grid::process_rank(i, j, k);
			if (face(4, process))
				ind_PWNEAZ = k;
			if (face(5, process))
				ind_PWFARZ = k;
		}

		if (i_less >= ind_PWNEAX + 1)
			nx[0] = fdtd::grid::process_rank(i_less, pos[1], pos[2]);
		else
			nx[0] = MPI_PROC_NULL;

		if (i_more <= ind_PWFARX - 1)
			nx[1] = fdtd::grid::process_rank(i_more, pos[1], pos[2]);
		else
			nx[1] = MPI_PROC_NULL;

		if (j_less >= ind_PWNEAY + 1)
			ny[0] = fdtd::grid::process_rank(pos[0], j_less, pos[2]);
		else
			ny[0] = MPI_PROC_NULL;

		if (j_more <= ind_PWFARY - 1)
			ny[1] = fdtd::grid::process_rank(pos[0], j_more, pos[2]);
		else
			ny[1] = MPI_PROC_NULL;

		if (k_less >= ind_PWNEAZ + 1)
			nz[0] = fdtd::grid::process_rank(pos[0], pos[1], k_less);
		else
			nz[0] = MPI_PROC_NULL;

		if (k_more <= ind_PWFARZ - 1)
			nz[1] = fdtd::grid::process_rank(pos[0], pos[1], k_more);
		else
			nz[1] = MPI_PROC_NULL;

		if (pos[0] <= ind_PWNEAX || pos[1] <= ind_PWNEAY || pos[2] <= ind_PWNEAZ || pos[0] >= ind_PWFARX || pos[1] >= ind_PWFARY || pos[2] >= ind_PWFARZ)
		{
			nx[0] = MPI_PROC_NULL;
			nx[1] = MPI_PROC_NULL;
			ny[0] = MPI_PROC_NULL;
			ny[1] = MPI_PROC_NULL;
			nz[0] = MPI_PROC_NULL;
			nz[1] = MPI_PROC_NULL;
		}
	}

	if (lat.periodic_x && lat.periodic_z)
	{

		i = int(0.5 * num_mpi_processes_x);
		k = int(0.5 * num_mpi_processes_z);

		for (j = 0; j < num_mpi_processes_y; ++j)
		{
			process = fdtd::grid::process_rank(i, j, k);
			if (face(2, process))
				ind_PWNEAY = j;
			if (face(3, process))
				ind_PWFARY = j;
		}

		if (i_less >= 0)
			nx[0] = fdtd::grid::process_rank(i_less, pos[1], pos[2]);
		else
		{
			if (lat.periodic_x)
				nx[0] = fdtd::grid::process_rank(num_mpi_processes_x - 1, pos[1], pos[2]);
			else
				nx[0] = MPI_PROC_NULL;
		}
		if (i_more < num_mpi_processes_x)
			nx[1] = fdtd::grid::process_rank(i_more, pos[1], pos[2]);
		else
		{
			if (lat.periodic_x)
				nx[1] = fdtd::grid::process_rank(0, pos[1], pos[2]);
			else
				nx[1] = MPI_PROC_NULL;
		}
		if (j_less >= ind_PWNEAY + 1)
			ny[0] = fdtd::grid::process_rank(pos[0], j_less, pos[2]);
		else
			ny[0] = MPI_PROC_NULL;

		if (j_more <= ind_PWFARY - 1)
			ny[1] = fdtd::grid::process_rank(pos[0], j_more, pos[2]);
		else
			ny[1] = MPI_PROC_NULL;

		if (k_less >= 0)
			nz[0] = fdtd::grid::process_rank(pos[0], pos[1], k_less);
		else
		{
			if (lat.periodic_z)
				nz[0] = fdtd::grid::process_rank(pos[0], pos[1], num_mpi_processes_z - 1);
			else
				nz[0] = MPI_PROC_NULL;
		}
		if (k_more < num_mpi_processes_z)
			nz[1] = fdtd::grid::process_rank(pos[0], pos[1], k_more);
		else
		{
			if (lat.periodic_z)
				nz[1] = fdtd::grid::process_rank(pos[0], pos[1], 0);
			else
				nz[1] = MPI_PROC_NULL;
		}

		if (pos[1] <= ind_PWNEAY || pos[1] >= ind_PWFARY)
		{
			nx[0] = MPI_PROC_NULL;
			nx[1] = MPI_PROC_NULL;
			ny[0] = MPI_PROC_NULL;
			ny[1] = MPI_PROC_NULL;
			nz[0] = MPI_PROC_NULL;
			nz[1] = MPI_PROC_NULL;
		}
	}

	return;
}

/**
 * @brief Uses the subdomain location indices to determine the rank of the
 * associated process
 * 
 * @param i x-index of the subdomain
 * @param j y-index of the subdomain
 * @param k z-index of the subdomain
 * @return int rank of process associated with this subdomain
 */
int fdtd::grid::process_rank(const int &i, const int &j, const int &k)
{
	static int N = num_mpi_processes_x * num_mpi_processes_y;

	if ((i >= 0 && i < num_mpi_processes_x) &&
		(j >= 0 && j < num_mpi_processes_y) &&
		(k >= 0 && k < num_mpi_processes_z))
	{

		return i + num_mpi_processes_x * j + N * k;
	}
	else
	{

		return -100;
	}
}

/**
 * @brief Get the Process Rank In Communicator object
 * 
 * @param tmpCommunicator 
 * @return int 
 */
int fdtd::grid::rank_in_communicator(MPI_Comm tmpCommunicator)
{
	int tmpRank;
	MPI_Comm_rank(tmpCommunicator, &tmpRank);
	return tmpRank;
}

/**
 * @brief Prints lenghts of subchunks along each axis and deletes temporary arrays to setup the decomposition
 * @return true when an error has occured
 * @return false when the end of the function is reached successfully
 */
bool fdtd::grid::print_decomposition(){
	//@gradient-decomposition
	auto list = [](const std::vector<int> &cells, int n) {
		std::string text;
		for (int i = 0; i < n; ++i)
			text += (i ? ", " : "") + std::to_string(cells[i]);
		return text;
	};
	if (fdtd::grid::world_rank == ROOTPROCESS)
	{
		fdtd::io::log("grid", "subdomain cells x: %s", list(fdtd::grid::cells_per_rank_x, num_mpi_processes_x).c_str());
		fdtd::io::log("grid", "subdomain cells y: %s", list(fdtd::grid::cells_per_rank_y, num_mpi_processes_y).c_str());
		fdtd::io::log("grid", "subdomain cells z: %s", list(fdtd::grid::cells_per_rank_z, num_mpi_processes_z).c_str());
	}
	//
	return false;
}
