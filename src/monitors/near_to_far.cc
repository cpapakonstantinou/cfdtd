/**
 * @file monitors/near_to_far.cc
 * @brief The near-to-far faces, monitors/near_to_far.h
 *
 * Replaces writeNearFar.cc: Save_Surface_Near_to_Far_YZ / XZ / XY, one copy per face orientation
 * and Write_Ini_File_NF.
 */
#include "sources/field_scale.h"

#include <complex>
#include <vector>
#include "monitors/near_to_far.h"
#include "monitors/state.h"
#include "grid/subdomain.h"
#include "core/views.h"
#include "grid/state.h"
#include "sources/state.h"
#include "io/h5io.h"
#include "monitors/flux_boxes.h"
#include <cstdio>
#include <string>
#include "grid/lattice.h"
#include "monitors/spec.h"
using fdtd::constants::C0;
using fdtd::constants::PI2;
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;

namespace
{
const fdtd::monitors::Spec &sp = fdtd::monitors::spec();        ///< what the monitors need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::monitors; // transitional: module state, monitors/state.h

namespace fdtd::monitors
{
namespace
{
/// the run parameters the transformation reads, binary, native layout
void write_parameters()
{
	const std::string path = sp.fields_dir + "/ini_near_to_far.dat";
	FILE *f = std::fopen(path.c_str(), "w+b");
	if (!f)
		return;
	const double cross_section = 1;
	const int n[3] = {lat.inner(0) - 2 * sp.nr_qsca[0], lat.inner(1) - 2 * sp.nr_qsca[1], lat.inner(2) - 2 * sp.nr_qsca[2]};
	const double center[3] = {double(sp.centre[0]), double(sp.centre[1] - sp.thickness / lat.d[1] / 2 - sp.antenna_shift / lat.d[1]),
							  double(sp.centre[2])};
	const int nea[3] = {sp.nr_qsca[0], sp.nr_qsca[1], sp.nr_qsca[2]};
	const int far[3] = {lat.inner(0) - sp.nr_qsca[0], lat.inner(1) - sp.nr_qsca[1], lat.inner(2) - sp.nr_qsca[2]};
	const double reference_wavelength = C0 / sp.fundamental_frequency;
	const double grid[7] = {lat.d[0], lat.d[1], lat.d[2], sp.dt, max_source_frequency, reference_wavelength, cross_section};
	std::fwrite(grid, sizeof(double), 7, f);
	std::fwrite(n, sizeof(int), 3, f);
	std::fwrite(center, sizeof(double), 3, f);
	std::fwrite(nea, sizeof(int), 3, f);
	std::fwrite(far, sizeof(int), 3, f);
	std::fwrite(&sp.polarization_deg, sizeof(int), 1, f);
	std::fclose(f);
}
} // namespace

void write_near_to_far()
{
	if (fdtd::grid::world_rank == ROOTPROCESS)
		write_parameters();

	static const char *axis_name[3] = {"x", "y", "z"};
	static const char *field_name[2][3] = {{"Ex", "Ey", "Ez"}, {"Hx", "Hy", "Hz"}};
	const int box_min[3] = {subdomain.i_qsca_min, subdomain.j_qsca_min, subdomain.k_qsca_min};
	const int box_max[3] = {subdomain.i_qsca_max, subdomain.j_qsca_max, subdomain.k_qsca_max};
	const bool on_face[3][2] = {{subdomain.qsca_nea_x, subdomain.qsca_far_x},
								{subdomain.qsca_nea_y, subdomain.qsca_far_y},
								{subdomain.qsca_nea_z, subdomain.qsca_far_z}};
	const MPI_Comm comm[3][2] = {{output_comm_sca_nea_x, output_comm_sca_far_x},
								 {output_comm_sca_nea_y, output_comm_sca_far_y},
								 {output_comm_sca_nea_z, output_comm_sca_far_z}};

	for (int n = 0; n < 3; ++n)
	{
		const int a = n == 0 ? 1 : 0;
		const int b = n == 2 ? 1 : 2; // the tangential axes, a < b
		if (!subdomain.cpml_only && box_min[a] != box_max[a] && box_min[b] != box_max[b])
			for (int m = 0; m < sp.nr_frequencies; ++m)
				for (int side = 0; side < 2; ++side)
				{
					if (!on_face[n][side])
						continue;
					const FluxFace &F = flux_boxes().box[0].face[n][side];
					const std::string path = sp.fields_dir + "/" + (side ? "far" : "near") + "_" + axis_name[n] + "_" +
											 std::to_string(std::lround((PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9))) + "nm.h5";
					const auto block = h5io::Block::plane(globalsizes_nf[a], globalsizes_nf[b], startindices_nf[a], startindices_nf[b],
														  localsize_nf[a], localsize_nf[b], subdomain.output_size[a], subdomain.output_size[b], box_min[a], box_min[b]);
					h5io::File file(path, comm[n][side]);
					file.attribute("domain_size", std::vector<int>{(fdtd::grid::lattice().inner(0) - 2 * sp.nr_qsca[0]), (fdtd::grid::lattice().inner(1) - 2 * sp.nr_qsca[1]), (fdtd::grid::lattice().inner(2) - 2 * sp.nr_qsca[2])});
					file.attribute("face", std::string(axis_name[n]));
					// the faces are stored normalised, normalize_plane; the harmonic band scaled as the field maps
					const fdtd::sources::FieldScale scale = fdtd::sources::field_scale(m, 0, phase_normalization, true);
					scale.note().write(file);
					std::vector<double> scaled_re;
					std::vector<double> scaled_im;
					for (int f = 0; f < 2; ++f)
						for (int c : {a, b})
						{
							const view4<double> re = F.re[f][0][c];
							const view4<double> im = F.im[f][0][c];
							const double *pr = &re[m][0][0][0];
							const double *pq = &im[m][0][0][0];
							if (scale.factor != std::complex<double>(1))
							{
								const std::size_t count = F.re[f][0][c].extent(1) * F.re[f][0][c].extent(2) * F.re[f][0][c].extent(3);
								scaled_re.resize(count);
								scaled_im.resize(count);
								for (std::size_t i = 0; i < count; ++i)
								{
									const std::complex<double> v = std::complex<double>(pr[i], pq[i]) * scale.factor;
									scaled_re[i] = v.real();
									scaled_im[i] = v.imag();
								}
								pr = scaled_re.data();
								pq = scaled_im.data();
							}
							h5io::write_complex(file, field_name[f][c], block, pr, pq);
						}
				}
		if (n < 2)
			MPI_Barrier(MPI_COMM_WORLD);
	}
}
} // namespace fdtd::monitors
