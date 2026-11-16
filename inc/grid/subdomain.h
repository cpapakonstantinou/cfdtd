#pragma once
/**
 * @file grid/subdomain.h
 * @brief One MPI rank's part of the lattice, module grid: this header, src/grid/subdomain.cc,
 *        the decomposition src/grid/decomposition.cc and the halo exchange src/grid/halo_exchange.cc.
 *
 * Every rank builds exactly one Subdomain, for itself, right after the decomposition:
 * its extents, its CPML ranges, the TF/SF and flux-box faces it owns, its plot planes.
 * Code that works on "this rank" reads `subdomain`, fdtd::grid::subdomain, a reference to this_rank(),
 * never a table indexed by rank.
 * fdtd::grid::set_lattice_data() fills it for this rank only, from the rank plan (rank_starts).
 */
#include <mpi.h>
#include <string>
#include "grid/planes.h"

namespace fdtd::grid
{
struct Subdomain
{
	int rank = -1;
	int position[3] = {};          ///< index in the rank grid along x, y, z
	int start_x = 0; ///< first global cell, PML included
	int start_y = 0; ///< first global cell, PML included
	int start_z = 0; ///< first global cell, PML included
	int end_x = 0; ///< one past the last global cell
	int end_y = 0; ///< one past the last global cell
	int end_z = 0; ///< one past the last global cell
	int cells_x = 0; ///< cells of this rank along each axis
	int cells_y = 0; ///< cells of this rank along each axis
	int cells_z = 0; ///< cells of this rank along each axis
	int owned_cells_x = 0; ///< cells updated here: one less on the last rank of an axis, shared far face
	int owned_cells_y = 0; ///< cells updated here: one less on the last rank of an axis, shared far face
	int owned_cells_z = 0; ///< cells updated here: one less on the last rank of an axis, shared far face
	double origin_x = 0; ///< [m] physical origin
	double origin_y = 0; ///< [m] physical origin
	double origin_z = 0; ///< [m] physical origin
	double length_x = 0; ///< [m] physical size
	double length_y = 0; ///< [m] physical size
	double length_z = 0; ///< [m] physical size
	int neighbours_x[2] = {}; ///< ranks before / after, MPI_PROC_NULL at a wall
	int neighbours_y[2] = {}; ///< ranks before / after, MPI_PROC_NULL at a wall
	int neighbours_z[2] = {}; ///< ranks before / after, MPI_PROC_NULL at a wall
	int neighbours_nl_x[2] = {}; ///< same, restricted to the TF/SF region, nonlinear fields
	int neighbours_nl_y[2] = {}; ///< same, restricted to the TF/SF region, nonlinear fields
	int neighbours_nl_z[2] = {}; ///< same, restricted to the TF/SF region, nonlinear fields
	bool cpml_only = {};
	int i_cpml_min = {};
	int j_cpml_min = {};
	int k_cpml_min = {};
	int i_cpml_max = {};
	int j_cpml_max = {};
	int k_cpml_max = {};
	// the CPML cells of this rank along each axis, in local indices: the near CPML is n <= *_nea_cpml_end,
	// the far CPML n >= *_far_cpml_start, which is also the origin of the far psi arrays, psi index
	// n - *_far_cpml_start; a rank without that CPML holds a value outside [0, cells]
	int i_nea_cpml_end = {};
	int j_nea_cpml_end = {};
	int k_nea_cpml_end = {};
	int i_far_cpml_start = {};
	int j_far_cpml_start = {};
	int k_far_cpml_start = {};
	int which_type_of_process_x = {};
	int which_type_of_process_y = {};
	int which_type_of_process_z = {};
	int shift_output_x = {};
	int shift_output_y = {};
	int shift_output_z = {};
	int inner_dim_x = {};
	int inner_dim_y = {};
	int inner_dim_z = {};
	int inner_org_x = {};
	int inner_org_y = {};
	int inner_org_z = {};
	int cpml_pad_nea_x = {};
	int cpml_pad_nea_y = {};
	int cpml_pad_nea_z = {};
	int cpml_pad_far_x = {};
	int cpml_pad_far_y = {};
	int cpml_pad_far_z = {};
	bool not_cpml_at_all = {};
	bool pw_nea_x = {};
	bool pw_far_x = {};
	bool pw_nea_y = {};
	bool pw_far_y = {};
	bool pw_nea_z = {};
	bool pw_far_z = {};
	int ipw_1 = {};
	int jpw_1 = {};
	int kpw_1 = {};
	int ipw_2 = {};
	int jpw_2 = {};
	int kpw_2 = {};
	int ipw_min = {};
	int ipw_max = {};
	int jpw_min = {};
	int jpw_max = {};
	int kpw_min = {};
	int kpw_max = {};
	int jpw_pad = {};
	int i_qsca_min = {};
	int i_qsca_max = {};
	int j_qsca_min = {};
	int j_qsca_max = {};
	int k_qsca_min = {};
	int k_qsca_max = {};
	int i_qtot_min = {};
	int i_qtot_max = {};
	int j_qtot_min = {};
	int j_qtot_max = {};
	int k_qtot_min = {};
	int k_qtot_max = {};
	bool qsca_nea_x = {};
	bool qsca_far_x = {};
	bool qsca_nea_y = {};
	bool qsca_far_y = {};
	bool qsca_nea_z = {};
	bool qsca_far_z = {};
	bool qtot_nea_x = {};
	bool qtot_far_x = {};
	bool qtot_nea_y = {};
	bool qtot_far_y = {};
	bool qtot_nea_z = {};
	bool qtot_far_z = {};
	int i_qsca_nea = {};
	int i_qsca_far = {};
	int j_qsca_nea = {};
	int j_qsca_far = {};
	int k_qsca_nea = {};
	int k_qsca_far = {};
	int i_qtot_nea = {};
	int i_qtot_far = {};
	int j_qtot_nea = {};
	int j_qtot_far = {};
	int k_qtot_nea = {};
	int k_qtot_far = {};
	bool animation_yz = {};
	int i_animation = {};
	bool art_nea_y = {};
	bool art_far_y = {};
	int j_art_1 = {};
	int j_art_2 = {};
	bool animation_xz = {};
	bool animation_xz_2 = {};
	bool animation_xz_1 = {};
	int j_animation = {};
	int j_animation_2 = {};
	int j_animation_1 = {};
	bool animation_xy = {};
	int k_animation = {};
	int shift_output_x_nf = {};
	int shift_output_y_nf = {};
	int shift_output_z_nf = {};
	int inner_org_x_nf = {};
	int inner_org_y_nf = {};
	int inner_org_z_nf = {};

	// set by app::initialise_run
	long output_start[3] = {};   ///< global index of this rank's first output cell, inner region, 0 on all-CPML ranks
	long output_size[3] = {};    ///< this rank's output cells per axis (inner_dim + 1), 0 on all-CPML ranks
	int j_for_uni = {};          ///< offset from this rank's j to the index on the 1D auxiliary grid
	bool has_plane_wave = {};    ///< this rank holds a TF/SF face of a plane-wave source

	// membership of the output communicators, was the ranks_excluded* lists
	bool nea_cpml_x() const { return which_type_of_process_x == 1 || which_type_of_process_x == 3 || which_type_of_process_x == 5; } ///< touches the near CPML along x
	bool far_cpml_x() const { return which_type_of_process_x == 2 || which_type_of_process_x == 4 || which_type_of_process_x == 5; } ///< touches the far CPML along x
	bool nea_cpml_y() const { return which_type_of_process_y == 1 || which_type_of_process_y == 3 || which_type_of_process_y == 5; } ///< touches the near CPML along y
	bool far_cpml_y() const { return which_type_of_process_y == 2 || which_type_of_process_y == 4 || which_type_of_process_y == 5; } ///< touches the far CPML along y
	bool nea_cpml_z() const { return which_type_of_process_z == 1 || which_type_of_process_z == 3 || which_type_of_process_z == 5; } ///< touches the near CPML along z
	bool far_cpml_z() const { return which_type_of_process_z == 2 || which_type_of_process_z == 4 || which_type_of_process_z == 5; } ///< touches the far CPML along z
	bool writes_output() const { return !cpml_only; }
	bool plots_yz() const { return !cpml_only && animation_yz; }
	bool plots_xz() const { return !cpml_only && animation_xz; }
	bool plots_xz_1() const { return !cpml_only && animation_xz_1; }
	bool plots_xz_2() const { return !cpml_only && animation_xz_2; }
	bool plots_xy() const { return !cpml_only && animation_xy; }
	bool flux_sca_nea_x() const { return !cpml_only && qsca_nea_x && j_qsca_min != j_qsca_max && k_qsca_min != k_qsca_max; }
	bool flux_sca_far_x() const { return !cpml_only && qsca_far_x && j_qsca_min != j_qsca_max && k_qsca_min != k_qsca_max; }
	bool flux_sca_nea_y() const { return !cpml_only && qsca_nea_y && i_qsca_min != i_qsca_max && k_qsca_min != k_qsca_max; }
	bool flux_sca_far_y() const { return !cpml_only && qsca_far_y && i_qsca_min != i_qsca_max && k_qsca_min != k_qsca_max; }
	bool flux_sca_nea_z() const { return !cpml_only && qsca_nea_z && i_qsca_min != i_qsca_max && j_qsca_min != j_qsca_max; }
	bool flux_sca_far_z() const { return !cpml_only && qsca_far_z && i_qsca_min != i_qsca_max && j_qsca_min != j_qsca_max; }
};

/// this rank's Subdomain, until the solver object owns it
Subdomain &this_rank();
inline Subdomain &subdomain = this_rank(); ///< this rank: the name the solver code uses

/// communicator of the ranks for which `member` is true, MPI_COMM_NULL elsewhere; MPI_COMM_WORLD when all are
MPI_Comm member_comm(bool member);
/// the lattice data of this rank: subdomain, CPML ranges, TF/SF, flux-box and plot-plane indices,
/// output communicators, decomposition.cc; true on failure, with why set to the reason
bool set_lattice_data(const GlobalPlanes &planes, std::string &why);
/// this rank's position (i, j, k) in the rank grid
void axis_indices(const int &rank, int pos[3]);
/// the neighbour ranks along x, y, z, -1 or MPI_PROC_NULL at the edges, periodic where set; _nl: for
/// the nonlocal margin exchange
void neighbour_ranks(int *pos, int *nx, int *ny, int *nz);
void neighbour_ranks_nl(int *pos, int *nx, int *ny, int *nz);
/// the rank at grid position (i, j, k)
int process_rank(const int &i, const int &j, const int &k);
/// this rank's rank in a communicator
int rank_in_communicator(MPI_Comm comm);
/// prints the cells per rank along each axis, root
bool print_decomposition();
} // namespace fdtd::grid
