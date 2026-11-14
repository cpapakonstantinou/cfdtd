/**
 * @file grid/state.cc
 * @brief Definitions of the grid module state, grid/state.h
 */
#include "grid/state.h"

namespace fdtd::grid
{
MPI_Datatype xy_interface_data_type_2 = MPI_DATATYPE_NULL;
MPI_Datatype xz_interface_data_type_2 = MPI_DATATYPE_NULL;
MPI_Datatype yz_interface_data_type_2 = MPI_DATATYPE_NULL;
int num_mpi_processes_x;
int num_mpi_processes_y;
int num_mpi_processes_z;
bool equal_decomp_x;
bool equal_decomp_y;
bool equal_decomp_z;
std::vector<int> arr_decomp_x;
std::vector<int> arr_decomp_y;
std::vector<int> arr_decomp_z;
std::vector<int> cells_per_rank_x;
std::vector<int> cells_per_rank_y;
std::vector<int> cells_per_rank_z;
int world_rank;
int world_size;
MPI_Datatype xy_interface_data_type = MPI_DATATYPE_NULL;
MPI_Datatype xz_interface_data_type = MPI_DATATYPE_NULL;
MPI_Datatype yz_interface_data_type = MPI_DATATYPE_NULL;
MPI_Comm output_comm;
MPI_Comm output_comm_plot_yz;
MPI_Comm output_comm_plot_xz;
MPI_Comm output_comm_plot_xz_2;
MPI_Comm output_comm_plot_xy;

} // namespace fdtd::grid
