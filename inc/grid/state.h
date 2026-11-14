/**
 * @file grid/state.h
 * @brief State of the grid module
 *
 * Moved out of the former shared namespace fdtd::state: read and written by this module, allocated or set by the
 * run set-up, app/setup.cc, and the configuration, applyConfig.cc. Defined in src/grid/state.cc.
 */
#ifndef FDTD_GRID_STATE_H
#define FDTD_GRID_STATE_H

#include "core/carray.h"
#include "core/constants.h"
#include <complex>
#include <mpi.h>
#include <string>
#include <vector>

namespace fdtd::grid
{
extern MPI_Datatype xy_interface_data_type_2;               ///< MPI datatype for one xy halo plane of a cells+3 array, exchange along z.
extern MPI_Datatype xz_interface_data_type_2;               ///< MPI datatype for one xz halo plane of a cells+3 array, exchange along y.
extern MPI_Datatype yz_interface_data_type_2;               ///< MPI datatype for one yz halo plane of a cells+3 array, exchange along x.
extern int num_mpi_processes_x;                             ///< Number of ranks along x.
extern int num_mpi_processes_y;                             ///< Number of ranks along y.
extern int num_mpi_processes_z;                             ///< Number of ranks along z.
extern bool equal_decomp_x;                                 ///< Split x evenly across ranks, otherwise use arr_decomp_x.
extern bool equal_decomp_y;                                 ///< Split y evenly across ranks, otherwise use arr_decomp_y.
extern bool equal_decomp_z;                                 ///< Split z evenly across ranks, otherwise use arr_decomp_z.
extern std::vector<int> arr_decomp_x;                       ///< Start offsets of the x slabs when the split is uneven.
extern std::vector<int> arr_decomp_y;                       ///< Start offsets of the y slabs when the split is uneven.
extern std::vector<int> arr_decomp_z;                       ///< Start offsets of the z slabs when the split is uneven.
extern std::vector<int> cells_per_rank_x;     ///< Interior cells of each rank along x, the decomposition.
extern std::vector<int> cells_per_rank_y;     ///< Interior cells of each rank along y, the decomposition.
extern std::vector<int> cells_per_rank_z;     ///< Interior cells of each rank along z, the decomposition.

extern int world_rank;                                       ///< Rank of this process in MPI_COMM_WORLD, set in main.
extern int world_size;                                       ///< Number of MPI processes, set in main.
extern MPI_Datatype xy_interface_data_type;                  ///< MPI datatype for one xy halo plane of H, exchange along z.
extern MPI_Datatype xz_interface_data_type;                  ///< MPI datatype for one xz halo plane of H, exchange along y.
extern MPI_Datatype yz_interface_data_type;                  ///< MPI datatype for one yz halo plane of H, exchange along x.
extern MPI_Comm output_comm;                                 ///< Communicator of the ranks that write volume output.
extern MPI_Comm output_comm_plot_yz;                         ///< Communicator of the ranks that cut the yz output plane.
extern MPI_Comm output_comm_plot_xz;                         ///< Communicator of the ranks that cut the xz output plane.
extern MPI_Comm output_comm_plot_xz_2;                       ///< Communicator of the ranks that cut xz output plane 2.
extern MPI_Comm output_comm_plot_xy;                         ///< Communicator of the ranks that cut the xy output plane.
} // namespace fdtd::grid

#endif
