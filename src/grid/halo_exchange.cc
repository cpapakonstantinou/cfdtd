/**
 * @file grid/halo_exchange.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-01-19
 * 
 * @brief Communicating across processes and initializing MPI structures
 * 
 * @details Houses functions for communicating fields over subdomains.
 * Ass well, it houses functions which create special custom MPI type
 * structures for probes and DFTs.
 * 
 * @copyright Copyright 2022
 */ 
#include "grid/subdomain.h"
#include "grid/lattice.h"
#include "grid/state.h"
#include "grid/halo_exchange.h"
// no shared state

namespace
{
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::grid; // transitional: module state, grid/state.h

namespace fdtd::grid
{
namespace
{
/// first element of plane idx normal to axis (0 x, 1 y, 2 z)
double *exchange_plane(view3<double> f, int axis, int idx) { return axis == 0 ? &f[idx][0][0] : axis == 1 ? &f[0][idx][0] : &f[0][0][idx]; }
} // namespace

/**
 * @brief Halo exchange of two field components along one axis
 * @verbatim
      E fields:  send plane hi to `more`, receive plane lo from `less`, the low halo plane is the
      H fields:  send plane lo to `less`, receive plane hi from `more`    neighbour's last plane
      lo = 0, hi = n                                   no periodic boundary along the axis
      lo = cpml, hi = n - cpml                         periodic, one rank along the axis
      periodic, several ranks: interior rank, less < more, lo = 0, hi = n; the first rank of the
      ring lo = cpml; the last one hi = n - cpml, the CPML slabs are skipped across the wrap
   @endverbatim
 * The first / last tests are the original rank-number tests, x: more % nx == 1 first, then
 * less % nx == nx - 2; z: less below the last layer first, then more beyond the first layer;
 * less and more may be the neighbours restricted to the TF/SF region, MPI_PROC_NULL at its edge.
 */
void exchange(Field field, int axis, view3<double> f1, view3<double> f2, int less, int more, int n, MPI_Datatype plane)
{
	const int tag = 0;
	MPI_Status status;
	const bool periodic = axis == 0 ? lat.periodic_x : axis == 2 ? lat.periodic_z : false;
	const int cpml = axis == 0 ? lat.cpml[0] : axis == 1 ? lat.cpml[1] : lat.cpml[2];
	const int ranks = axis == 0 ? num_mpi_processes_x : axis == 1 ? num_mpi_processes_y : num_mpi_processes_z;
	int lo = 0, hi = n;
	if (periodic)
	{
		if (ranks == 1)
			lo = cpml, hi = n - cpml;
		else if (less >= more)
		{
			const int layer = num_mpi_processes_x * num_mpi_processes_y;
			const bool first = axis == 0 ? more % num_mpi_processes_x == 1 : !(less < fdtd::grid::world_size - layer) && more >= layer;
			const bool last = axis == 0 ? !first && less % num_mpi_processes_x == num_mpi_processes_x - 2 : less < fdtd::grid::world_size - layer;
			if (first)
				lo = cpml;
			else if (last)
				hi = n - cpml;
			else
				return; // as before: no exchange
		}
	}
	for (view3<double> f : {f1, f2})
	{
		if (field == Field::E)
			MPI_Sendrecv(exchange_plane(f, axis, hi), 1, plane, more, tag, exchange_plane(f, axis, lo), 1, plane, less, tag, MPI_COMM_WORLD, &status);
		else
			MPI_Sendrecv(exchange_plane(f, axis, lo), 1, plane, less, tag, exchange_plane(f, axis, hi), 1, plane, more, tag, MPI_COMM_WORLD, &status);
	}
}
} // namespace fdtd::grid

/**
 * @brief communicating several boundary vectors between processes
 * 
 * @param VecX X vector to be shared
 * @param VecY Y vector to be shared
 * @param VecZ Z vector to be shared
 * 
 * @details 
 */
void fdtd::grid::exchange_margins(view3<double> VecX, view3<double> VecY, view3<double> VecZ)
{
	int NX, NY, NZ;
	int tag = 0;
	MPI_Status status;
	int moreX, lessX;
	int moreZ, lessZ;
	int moreY, lessY;
	NX = subdomain.cells_x;
	NY = subdomain.cells_y;
	NZ = subdomain.cells_z;
	lessX = subdomain.neighbours_x[0];
	moreX = subdomain.neighbours_x[1];
	lessY = subdomain.neighbours_y[0];
	moreY = subdomain.neighbours_y[1];
	lessZ = subdomain.neighbours_z[0];
	moreZ = subdomain.neighbours_z[1];

	if (lat.periodic_x)
	{
		if (num_mpi_processes_x == 1)
		{
			MPI_Sendrecv(&VecX[NX - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag,
									 &VecX[0 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status); // NX-1 -> -1
			MPI_Sendrecv(&VecX[1 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag,
									 &VecX[NX + 1 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status); // 0 and 1 -> NX and NX+1

			MPI_Sendrecv(&VecX[2 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag,
									 &VecX[NX + 2 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
		}

		else
		{
			if (lessX < moreX)
			{
				MPI_Sendrecv(&VecX[NX][0][0], 1, yz_interface_data_type_2, moreX, tag,
										 &VecX[0][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
				MPI_Sendrecv(&VecX[1][0][0], 1, yz_interface_data_type_2, lessX, tag,
										 &VecX[NX + 1][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
				MPI_Sendrecv(&VecX[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
										 &VecX[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
			}
			else
			{
				if (moreX % num_mpi_processes_x == 1)
				{
					MPI_Sendrecv(&VecX[NX][0][0], 1, yz_interface_data_type_2, moreX, tag,
											 &VecX[0 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecX[1 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag,
											 &VecX[NX + 1][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecX[2 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag,
											 &VecX[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
				}
				else if (lessX % num_mpi_processes_x == num_mpi_processes_x - 2)

				{
					MPI_Sendrecv(&VecX[NX - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag,
											 &VecX[0][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecX[1][0][0], 1, yz_interface_data_type_2, lessX, tag,
											 &VecX[NX + 1 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecX[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
											 &VecX[NX + 2 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
				}
			}
		}
	}
	else
	{
		MPI_Sendrecv(&VecX[NX][0][0], 1, yz_interface_data_type_2, moreX, tag,
								 &VecX[0][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&VecX[1][0][0], 1, yz_interface_data_type_2, lessX, tag,
								 &VecX[NX + 1][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&VecX[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
								 &VecX[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
	}

	//Y-Interface
	// MPI_Sendrecv(&px_nl_2[0][NY][0], 1, XZ_INTERFACE_DATA_TYPE_3, moreY, tag, // NY-1 and NY -> -1 and 0
	MPI_Sendrecv(&VecX[0][2][0], 1, xz_interface_data_type_2, lessY, tag, // 1 -> NY+1
							 &VecX[0][NY + 2][0], 1, xz_interface_data_type_2, moreY, tag, MPI_COMM_WORLD, &status);

	MPI_Sendrecv(&VecX[0][NY + 1][0], 1, xz_interface_data_type_2, moreY, tag, // NY -> 0
							 &VecX[0][1][0], 1, xz_interface_data_type_2, lessY, tag, MPI_COMM_WORLD, &status);

	//Z-Interface
	if (lat.periodic_z)
	{
		if (num_mpi_processes_z == 1)
		{
			// MPI_Sendrecv(&px_nl_2[0][0][NZ-num_cpml_cells_z], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag, //NZ-1 and NZ -> -1 and 0
			MPI_Sendrecv(&VecX[0][0][2 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag, // 1 -> NZ+1
									 &VecX[0][0][NZ + 2 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(&VecX[0][0][NZ + 1 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, //NZ -> 0
									 &VecX[0][0][1 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
		}
		else
		{
			if (lessZ < moreZ)
			{
				// MPI_Sendrecv(&px_nl_2[0][0][NZ], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
				MPI_Sendrecv(&VecX[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
										 &VecX[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);

				MPI_Sendrecv(&VecX[0][0][NZ + 1], 1, xy_interface_data_type_2, moreZ, tag,
										 &VecX[0][0][1], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
			}
			else
			{
				if (lessZ < fdtd::grid::world_size - num_mpi_processes_x * num_mpi_processes_y)
				{
					// MPI_Sendrecv(&px_nl_2[0][0][NZ-num_cpml_cells_z], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
					MPI_Sendrecv(&VecX[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
											 &VecX[0][0][NZ + 2 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecX[0][0][NZ + 1 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag,
											 &VecX[0][0][1], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
				}
				else if (moreZ >= num_mpi_processes_x * num_mpi_processes_y)
				{
					// MPI_Sendrecv(&px_nl_2[0][0][NZ], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
					MPI_Sendrecv(&VecX[0][0][2 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag,
											 &VecX[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);

					MPI_Sendrecv(&VecX[0][0][NZ + 1], 1, xy_interface_data_type_2, moreZ, tag,
											 &VecX[0][0][1 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
				}
			}
		}
	}
	else
	{
		// MPI_Sendrecv(&px_nl_2[0][0][NZ], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
		MPI_Sendrecv(&VecX[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
								 &VecX[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&VecX[0][0][NZ + 1], 1, xy_interface_data_type_2, moreZ, tag,
								 &VecX[0][0][1], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////

	//Py:
	//X-interface
	if (lat.periodic_x)
	{
		if (num_mpi_processes_x == 1)
		{
			// MPI_Sendrecv(&py_nl_2[NX-num_cpml_cells_x][0][0], 1, YZ_INTERFACE_DATA_TYPE_3, moreX, tag, // Nx-1 and Nx -> -1 and 0
			MPI_Sendrecv(&VecY[2 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag, // 1 -> Nx+1
									 &VecY[NX + 2 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(&VecY[NX + 1 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, // Nx -> 0
									 &VecY[1 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
		}

		else
		{
			if (lessX < moreX)
			{
				// 	MPI_Sendrecv(&py_nl_2[NX][0][0], 1, YZ_INTERFACE_DATA_TYPE_3, moreX, tag,
				MPI_Sendrecv(&VecY[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
										 &VecY[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
				MPI_Sendrecv(&VecY[NX + 1][0][0], 1, yz_interface_data_type_2, moreX, tag,
										 &VecY[1][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
			}
			else
			{
				if (moreX % num_mpi_processes_x == 1)
				{
					// MPI_Sendrecv(&py_nl_2[NX][0][0], 1, YZ_INTERFACE_DATA_TYPE_3, moreX, tag,
					MPI_Sendrecv(&VecY[2 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag,
											 &VecY[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecY[NX + 1][0][0], 1, yz_interface_data_type_2, moreX, tag,
											 &VecY[1 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
				}
				else if (lessX % num_mpi_processes_x == num_mpi_processes_x - 2)

				{
					// MPI_Sendrecv(&py_nl_2[NX-num_cpml_cells_x][0][0], 1, YZ_INTERFACE_DATA_TYPE_3, moreX, tag,
					MPI_Sendrecv(&VecY[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
											 &VecY[NX + 2 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecY[NX + 1 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag,
											 &VecY[1][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
				}
			}
		}
	}
	else
	{
		// MPI_Sendrecv(&py_nl_2[NX][0][0], 1, YZ_INTERFACE_DATA_TYPE_3, moreX, tag,
		MPI_Sendrecv(&VecY[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
								 &VecY[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&VecY[NX + 1][0][0], 1, yz_interface_data_type_2, moreX, tag,
								 &VecY[1][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
	}

	//Y-Interface
	// MPI_Sendrecv(&py_nl_2[0][NY][0], 1, xz_interface_data_type_2, moreY, tag, // NY-1 -> -1
	// MPI_Sendrecv(&py_nl_2[0][1][0], 1, XZ_INTERFACE_DATA_TYPE_3, lessY, tag,          //0 and 1 -> NY and NY+1

	MPI_Sendrecv(&VecY[0][NY][0], 1, xz_interface_data_type_2, moreY, tag, // NY-1 -> -1
							 &VecY[0][0][0], 1, xz_interface_data_type_2, lessY, tag, MPI_COMM_WORLD, &status);
	MPI_Sendrecv(&VecY[0][1][0], 1, xz_interface_data_type_2, lessY, tag, //0  -> NY
							 &VecY[0][NY + 1][0], 1, xz_interface_data_type_2, moreY, tag, MPI_COMM_WORLD, &status);
	MPI_Sendrecv(&VecY[0][2][0], 1, xz_interface_data_type_2, lessY, tag, //1  -> NY + 1
							 &VecY[0][NY + 2][0], 1, xz_interface_data_type_2, moreY, tag, MPI_COMM_WORLD, &status);

	//Z-Interface
	if (lat.periodic_z)
	{
		if (num_mpi_processes_z == 1)
		{
			// MPI_Sendrecv(&py_nl_2[0][0][NZ-num_cpml_cells_z], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag, //NZ-1 and NZ -> -1 and 0
			MPI_Sendrecv(&VecY[0][0][2 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag, //1 -> NZ+1
									 &VecY[0][0][NZ + 2 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);

			MPI_Sendrecv(&VecY[0][0][NZ + 1 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, // NZ ->  0
									 &VecY[0][0][1 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
		}
		else
		{
			if (lessZ < moreZ)
			{
				// MPI_Sendrecv(&py_nl_2[0][0][NZ], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
				MPI_Sendrecv(&VecY[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
										 &VecY[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);

				MPI_Sendrecv(&VecY[0][0][NZ + 1], 1, xy_interface_data_type_2, moreZ, tag,
										 &VecY[0][0][1], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
			}
			else
			{
				if (lessZ < fdtd::grid::world_size - num_mpi_processes_x * num_mpi_processes_y)
				{
					// MPI_Sendrecv(&py_nl_2[0][0][NZ-num_cpml_cells_z], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
					MPI_Sendrecv(&VecY[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
											 &VecY[0][0][NZ + 2 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecY[0][0][NZ + 1 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag,
											 &VecY[0][0][1], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
				}
				else if (moreZ >= num_mpi_processes_x * num_mpi_processes_y)
				{
					// MPI_Sendrecv(&py_nl_2[0][0][NZ], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
					MPI_Sendrecv(&VecY[0][0][2 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag,
											 &VecY[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecY[0][0][NZ + 1], 1, xy_interface_data_type_2, moreZ, tag,
											 &VecY[0][0][1 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
				}
			}
		}
	}
	else
	{
		// MPI_Sendrecv(&py_nl_2[0][0][NZ], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
		MPI_Sendrecv(&VecY[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
								 &VecY[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&VecY[0][0][NZ + 1], 1, xy_interface_data_type_2, moreZ, tag,
								 &VecY[0][0][1], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
	//Pz
	//X-interface

	if (lat.periodic_x)
	{
		if (num_mpi_processes_x == 1)
		{
			// MPI_Sendrecv(&pz_nl_2[NX-num_cpml_cells_x][0][0], 1, YZ_INTERFACE_DATA_TYPE_3, moreX, tag,  // NX-1 and NX -> -1 and 0
			MPI_Sendrecv(&VecZ[2 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag,
									 &VecZ[NX + 2 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(&VecZ[NX + 1 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, // NX ->  0
									 &VecZ[1 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
		}

		else
		{
			if (lessX < moreX)
			{
				// 	MPI_Sendrecv(&pz_nl_2[NX][0][0], 1, YZ_INTERFACE_DATA_TYPE_3, moreX, tag,
				MPI_Sendrecv(&VecZ[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
										 &VecZ[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
				MPI_Sendrecv(&VecZ[NX + 1][0][0], 1, yz_interface_data_type_2, moreX, tag,
										 &VecZ[1][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
			}
			else
			{
				if (moreX % num_mpi_processes_x == 1)
				{
					// MPI_Sendrecv(&pz_nl_2[NX][0][0], 1, YZ_INTERFACE_DATA_TYPE_3, moreX, tag,
					MPI_Sendrecv(&VecZ[2 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag,
											 &VecZ[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecZ[NX + 1][0][0], 1, yz_interface_data_type_2, moreX, tag,
											 &VecZ[1 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
				}
				else if (lessX % num_mpi_processes_x == num_mpi_processes_x - 2)

				{
					// MPI_Sendrecv(&pz_nl_2[NX-num_cpml_cells_x][0][0], 1, YZ_INTERFACE_DATA_TYPE_3, moreX, tag,
					MPI_Sendrecv(&VecZ[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
											 &VecZ[NX + 2 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecZ[NX + 1 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag,
											 &VecZ[1][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
				}
			}
		}
	}
	else
	{
		// MPI_Sendrecv(&pz_nl_2[NX][0][0], 1, YZ_INTERFACE_DATA_TYPE_3, moreX, tag,
		MPI_Sendrecv(&VecZ[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
								 &VecZ[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&VecZ[NX + 1][0][0], 1, yz_interface_data_type_2, moreX, tag,
								 &VecZ[1][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
	}

	//Y-Interface
	// MPI_Sendrecv(&pz_nl_2[0][NY][0], 1, XZ_INTERFACE_DATA_TYPE_3, moreY, tag, //NY-1 and NY -> -1 and 0
	MPI_Sendrecv(&VecZ[0][2][0], 1, xz_interface_data_type_2, lessY, tag,
							 &VecZ[0][NY + 2][0], 1, xz_interface_data_type_2, moreY, tag, MPI_COMM_WORLD, &status);
	MPI_Sendrecv(&VecZ[0][NY + 1][0], 1, xz_interface_data_type_2, moreY, tag, //NY -> 0
							 &VecZ[0][1][0], 1, xz_interface_data_type_2, lessY, tag, MPI_COMM_WORLD, &status);
	//Z-Interface

	if (lat.periodic_z)
	{
		if (num_mpi_processes_z == 1)
		{
			// MPI_Sendrecv(&pz_nl_2[0][0][NZ-num_cpml_cells_z], 1, xy_interface_data_type_2, moreZ, tag,
			// MPI_Sendrecv(&pz_nl_2[0][0][1+num_cpml_cells_z], 1, XY_INTERFACE_DATA_TYPE_3, lessZ, tag,

			MPI_Sendrecv(&VecZ[0][0][NZ - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag,
									 &VecZ[0][0][0 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(&VecZ[0][0][1 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag,
									 &VecZ[0][0][NZ + 1 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(&VecZ[0][0][2 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag,
									 &VecZ[0][0][NZ + 2 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
		}
		else
		{
			if (lessZ < moreZ)
			{
				// MPI_Sendrecv(&pz_nl_2[0][0][NZ], 1, xy_interface_data_type_2, moreZ, tag,
				// MPI_Sendrecv(&pz_nl_2[0][0][1], 1, XY_INTERFACE_DATA_TYPE_3, lessZ, tag,

				MPI_Sendrecv(&VecZ[0][0][NZ], 1, xy_interface_data_type_2, moreZ, tag,
										 &VecZ[0][0][0], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status); //NZ-1 -> -1
				MPI_Sendrecv(&VecZ[0][0][1], 1, xy_interface_data_type_2, lessZ, tag,
										 &VecZ[0][0][NZ + 1], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status); //0  -> NZ
				MPI_Sendrecv(&VecZ[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
										 &VecZ[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status); //1  -> NZ+1
			}
			else
			{
				if (lessZ < fdtd::grid::world_size - num_mpi_processes_x * num_mpi_processes_y)
				{
					// MPI_Sendrecv(&pz_nl_2[0][0][NZ-num_cpml_cells_z], 1, xy_interface_data_type_2, moreZ, tag,
					// MPI_Sendrecv(&pz_nl_2[0][0][1], 1, XY_INTERFACE_DATA_TYPE_3, lessZ, tag,

					MPI_Sendrecv(&VecZ[0][0][NZ - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag,
											 &VecZ[0][0][0], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecZ[0][0][1], 1, xy_interface_data_type_2, lessZ, tag,
											 &VecZ[0][0][NZ + 1 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecZ[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
											 &VecZ[0][0][NZ + 2 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
				}
				else if (moreZ >= num_mpi_processes_x * num_mpi_processes_y)
				{
					// MPI_Sendrecv(&pz_nl_2[0][0][NZ], 1, xy_interface_data_type_2, moreZ, tag,
					// MPI_Sendrecv(&pz_nl_2[0][0][1+num_cpml_cells_z], 1, XY_INTERFACE_DATA_TYPE_3, lessZ, tag,
					MPI_Sendrecv(&VecZ[0][0][NZ], 1, xy_interface_data_type_2, moreZ, tag,
											 &VecZ[0][0][0 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecZ[0][0][1 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag,
											 &VecZ[0][0][NZ + 1], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecZ[0][0][2 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag,
											 &VecZ[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
				}
			}
		}
	}
	else
	{
		// MPI_Sendrecv(&pz_nl_2[0][0][NZ], 1, xy_interface_data_type_2, moreZ, tag,
		// MPI_Sendrecv(&pz_nl_2[0][0][1], 1, XY_INTERFACE_DATA_TYPE_3, lessZ, tag,
		MPI_Sendrecv(&VecZ[0][0][NZ], 1, xy_interface_data_type_2, moreZ, tag,
								 &VecZ[0][0][0], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&VecZ[0][0][1], 1, xy_interface_data_type_2, lessZ, tag,
								 &VecZ[0][0][NZ + 1], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&VecZ[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
								 &VecZ[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
	}
}

/**
 * @brief 
 * 
 * @param VecX 
 */
void fdtd::grid::exchange_margin(view3<double> VecX)
{
	int NX, NY, NZ;
	int tag = 0;
	MPI_Status status;
	int moreX, lessX;
	int moreZ, lessZ;
	int moreY, lessY;
	NX = subdomain.cells_x;
	NY = subdomain.cells_y;
	NZ = subdomain.cells_z;
	lessX = subdomain.neighbours_x[0];
	moreX = subdomain.neighbours_x[1];
	lessY = subdomain.neighbours_y[0];
	moreY = subdomain.neighbours_y[1];
	lessZ = subdomain.neighbours_z[0];
	moreZ = subdomain.neighbours_z[1];

	if (lat.periodic_x)
	{
		if (num_mpi_processes_x == 1)
		{
			MPI_Sendrecv(&VecX[NX - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag,
									 &VecX[0 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status); // NX-1 -> -1

			MPI_Sendrecv(&VecX[2 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag,
									 &VecX[NX + 2 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status); // 0 and 1 -> NX and NX+1
		}

		else
		{
			if (lessX < moreX)
			{
				MPI_Sendrecv(&VecX[NX][0][0], 1, yz_interface_data_type_2, moreX, tag,
										 &VecX[0][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
				MPI_Sendrecv(&VecX[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
										 &VecX[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
			}
			else
			{
				if (moreX % num_mpi_processes_x == 1)
				{
					MPI_Sendrecv(&VecX[NX][0][0], 1, yz_interface_data_type_2, moreX, tag,
											 &VecX[0 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecX[2 + lat.cpml[0]][0][0], 1, yz_interface_data_type_2, lessX, tag,
											 &VecX[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
				}
				else if (lessX % num_mpi_processes_x == num_mpi_processes_x - 2)

				{
					MPI_Sendrecv(&VecX[NX - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag,
											 &VecX[0][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecX[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
											 &VecX[NX + 2 - lat.cpml[0]][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
				}
			}
		}
	}
	else
	{
		MPI_Sendrecv(&VecX[NX][0][0], 1, yz_interface_data_type_2, moreX, tag,
								 &VecX[0][0][0], 1, yz_interface_data_type_2, lessX, tag, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&VecX[2][0][0], 1, yz_interface_data_type_2, lessX, tag,
								 &VecX[NX + 2][0][0], 1, yz_interface_data_type_2, moreX, tag, MPI_COMM_WORLD, &status);
	}

	//Y-Interface
	// MPI_Sendrecv(&px_nl_2[0][NY][0], 1, XZ_INTERFACE_DATA_TYPE_3, moreY, tag, // NY-1 and NY -> -1 and 0
	// MPI_Sendrecv(&px_nl_2[0][2][0], 1, xz_interface_data_type_2, lessY, tag, // 1 -> NY+1

	MPI_Sendrecv(&VecX[0][NY][0], 1, xz_interface_data_type_2, moreY, tag, // NY -> 0
							 &VecX[0][0][0], 1, xz_interface_data_type_2, lessY, tag, MPI_COMM_WORLD, &status);

	MPI_Sendrecv(&VecX[0][2][0], 1, xz_interface_data_type_2, lessY, tag, // NY -> 0
							 &VecX[0][NY + 2][0], 1, xz_interface_data_type_2, moreY, tag, MPI_COMM_WORLD, &status);

	//Z-Interface
	if (lat.periodic_z)
	{
		if (num_mpi_processes_z == 1)
		{
			// MPI_Sendrecv(&px_nl_2[0][0][NZ-num_cpml_cells_z], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag, //NZ-1 and NZ -> -1 and 0

			MPI_Sendrecv(&VecX[0][0][2 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag, // 1 -> NZ+1
									 &VecX[0][0][NZ + 2 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(&VecX[0][0][NZ - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, //NZ -> 0
									 &VecX[0][0][0 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
		}
		else
		{
			if (lessZ < moreZ)
			{
				// MPI_Sendrecv(&px_nl_2[0][0][NZ], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
				MPI_Sendrecv(&VecX[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
										 &VecX[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);

				MPI_Sendrecv(&VecX[0][0][NZ], 1, xy_interface_data_type_2, moreZ, tag,
										 &VecX[0][0][0], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
			}
			else
			{
				if (lessZ < fdtd::grid::world_size - num_mpi_processes_x * num_mpi_processes_y)
				{
					// MPI_Sendrecv(&px_nl_2[0][0][NZ-num_cpml_cells_z], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
					MPI_Sendrecv(&VecX[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
											 &VecX[0][0][NZ + 2 - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecX[0][0][NZ - lat.cpml[2]], 1, xy_interface_data_type_2, moreZ, tag,
											 &VecX[0][0][0], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
				}
				else if (moreZ >= num_mpi_processes_x * num_mpi_processes_y)
				{
					// MPI_Sendrecv(&px_nl_2[0][0][NZ], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
					MPI_Sendrecv(&VecX[0][0][2 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag,
											 &VecX[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
					MPI_Sendrecv(&VecX[0][0][NZ], 1, xy_interface_data_type_2, moreZ, tag,
											 &VecX[0][0][0 + lat.cpml[2]], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
				}
			}
		}
	}
	else
	{
		// MPI_Sendrecv(&px_nl_2[0][0][NZ], 1, XY_INTERFACE_DATA_TYPE_3, moreZ, tag,
		MPI_Sendrecv(&VecX[0][0][2], 1, xy_interface_data_type_2, lessZ, tag,
								 &VecX[0][0][NZ + 2], 1, xy_interface_data_type_2, moreZ, tag, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&VecX[0][0][NZ], 1, xy_interface_data_type_2, moreZ, tag,
								 &VecX[0][0][0], 1, xy_interface_data_type_2, lessZ, tag, MPI_COMM_WORLD, &status);
	}
}

