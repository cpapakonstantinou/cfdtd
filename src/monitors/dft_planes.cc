/**
 * @file monitors/dft_planes.cc
 * @brief The display DFT planes, monitors/dft_planes.h
 */
#include "monitors/dft_planes.h"
#include "monitors/state.h"
#include "grid/subdomain.h"
#include "core/views.h"
#include "grid/state.h"
#include "sources/state.h"
#include "monitors/plane_dft.h"
#include "io/h5io.h"
#include <cmath>
#include <string>
#include <vector>
#include "monitors/spec.h"
#include "grid/lattice.h"
using fdtd::constants::C0;
using fdtd::constants::PI2;
using fdtd::grid::subdomain;

namespace
{
const fdtd::monitors::Spec &sp = fdtd::monitors::spec();        ///< what the monitors need from the run
} // namespace
using namespace fdtd::monitors; // transitional: module state, monitors/state.h

namespace fdtd::monitors
{
DftPlanes &dft_planes()
{
	static DftPlanes p;
	return p;
}

namespace
{
bool plane_on(int n)
{
	return n == 0 ? subdomain.plots_yz() && dft_plane_yz : n == 1 ? subdomain.plots_xz() && (dft_plane_xz || sp.graphene_optim) : subdomain.plots_xy() && dft_plane_xy;
}

/// allocates fields f (0 E, 1 H) of set s of a plane with normal n
void allocate_plane(DftPlane &P, int n, int f, int s)
{
	const int dims[3] = {n == 0 ? 0 : subdomain.inner_dim_x, n == 1 ? 0 : subdomain.inner_dim_y, n == 2 ? 0 : subdomain.inner_dim_z};
	for (int c = 0; c < 3; ++c)
	{
		P.re[f][s][c] = field4<double>(sp.nr_frequencies, dims[0], dims[1], dims[2]);
		P.im[f][s][c] = field4<double>(sp.nr_frequencies, dims[0], dims[1], dims[2]);
	}
}

void accumulate_plane(DftPlane &P, int f, int s, int n, double t, int p)
{
	sample_plane(n, p, f, s, t, P.re[f][s], P.im[f][s]);
}
} // namespace

void DftPlanes::allocate()
{
	if (subdomain.cpml_only)
		return;
	for (int n = 0; n < 3; ++n)
		if (plane_on(n))
			for (int s = 0; s < (sp.hydrodynamics ? 2 : 1); ++s)
				for (int f = 0; f < 2; ++f)
					allocate_plane(plane[n], n, f, s);
	if (subdomain.plots_xz() && (dft_plane_xz || sp.graphene_optim) && sp.graphene_optim)
		allocate_plane(graphene, 1, 0, 0);
	if (subdomain.plots_xz_1() && dft_xz_1)
		allocate_plane(xz_1, 1, 0, 0);
	if (subdomain.plots_xz_2() && dft_xz_2)
		allocate_plane(xz_2, 1, 0, 0);
}

void DftPlanes::release()
{
	for (DftPlane *P : {&plane[0], &plane[1], &plane[2], &xz_1, &xz_2, &graphene})
		for (int f = 0; f < 2; ++f)
			for (int s = 0; s < 2; ++s)
				for (int c = 0; c < 3; ++c)
				{
					P->re[f][s][c].reset();
					P->im[f][s][c].reset();
				}
}

/**
 * Order as before: the xy plane, E, H, harmonic E, harmonic H; the xz plane, E, H, the graphene
 * plane, harmonic E, harmonic H; the forward and backward far-field planes; the yz plane.
 */
void DftPlanes::sample(double t)
{
	if (subdomain.cpml_only)
		return;
	auto full = [&](int n, int p, bool with_graphene) {
		accumulate_plane(plane[n], 0, 0, n, t, p);
		accumulate_plane(plane[n], 1, 0, n, t, p);
		if (with_graphene)
			accumulate_plane(graphene, 0, 0, 1, t, p - 1);
		if (sp.hydrodynamics)
		{
			accumulate_plane(plane[n], 0, 1, n, t, p);
			accumulate_plane(plane[n], 1, 1, n, t, p);
		}
	};
	if (plane_on(2))
		full(2, subdomain.k_animation, false);
	if (plane_on(1))
		full(1, subdomain.j_animation, sp.graphene_optim);
	if (subdomain.plots_xz_2() && dft_xz_2)
		accumulate_plane(xz_2, 0, 0, 1, t, subdomain.j_animation_2);
	if (subdomain.plots_xz_1() && dft_xz_1)
		accumulate_plane(xz_1, 0, 0, 1, t, subdomain.j_animation_1);
	if (plane_on(0))
		full(0, subdomain.i_animation, false);
}

void DftPlanes::normalize()
{
	if (subdomain.cpml_only)
		return;
	for (int n = 0; n < 3; ++n)
		if (plane_on(n))
			for (int f = 0; f < 2; ++f)
				normalize_plane(n, plane[n].re[f][0], plane[n].im[f][0]);
	if (plane_on(1) && sp.graphene_optim)
		normalize_plane(1, graphene.re[0][0], graphene.im[0][0]);
	if (subdomain.plots_xz_2() && dft_xz_2)
		normalize_plane(1, xz_2.re[0][0], xz_2.im[0][0]);
	if (subdomain.plots_xz_1() && dft_xz_1)
		normalize_plane(1, xz_1.re[0][0], xz_1.im[0][0]);
}

namespace
{
/**
 * One plane DFT at frequency m into an HDF5 file. Each rank holds a, subdomain.output_size[a] x subdomain.output_size[b]
 * buffer whose last row/column is a ghost cell; the, subdomain.output_size - 1, interior is placed at
 * subdomain.output_start in a global (globalsizes[a] - 1) x (globalsizes[b] - 1) plane.
 */
void write_map(const std::string &path, MPI_Comm comm, int a, int b, int m, const ctetrad<double> (&re)[3], const ctetrad<double> (&im)[3])
{
	const auto block = h5io::Block::plane(sp.global_size[a] - 1, sp.global_size[b] - 1, subdomain.output_start[a], subdomain.output_start[b], subdomain.output_size[a] - 1,
										  subdomain.output_size[b] - 1, subdomain.output_size[a], subdomain.output_size[b]);
	h5io::File file(path, comm);
	file.attribute("domain_size", std::vector<int>{fdtd::grid::lattice().inner(0), fdtd::grid::lattice().inner(1), fdtd::grid::lattice().inner(2)});
	file.attribute("plane", std::string{"XYZ"[a], "XYZ"[b]});
	file.attribute("wavelength", double(PI2 * C0 / fdtd::sources::analysis_omega[m]));
	for (int c = 0; c < 3; ++c)
	{
		const view4<double> r = re[c], q = im[c];
		h5io::write_complex(file, std::string(1, "xyz"[c]), block, &r[m][0][0][0], &q[m][0][0][0]);
	}
}
} // namespace

void DftPlanes::write_maps()
{
	if (subdomain.cpml_only)
		return;
	const std::string dir = sp.fields_dir + "/";
	static const char *tag[3] = {"YZ", "XZ", "XY"};
	const MPI_Comm comm[3] = {fdtd::grid::output_comm_plot_yz, fdtd::grid::output_comm_plot_xz, fdtd::grid::output_comm_plot_xy};
	const int axes[3][2] = {{1, 2}, {0, 2}, {0, 1}};
	// E, H and their harmonic sets of the plane with normal n, as before: E, H, E_nonlinear, H_nonlinear
	auto full = [&](int n, int m, const std::string &wl) {
		const std::string suffix = wl + "_" + tag[n] + ".h5";
		const int a = axes[n][0], b = axes[n][1];
		write_map(dir + "E_" + suffix, comm[n], a, b, m, plane[n].re[0][0], plane[n].im[0][0]);
		if (dft_h)
			write_map(dir + "H_" + suffix, comm[n], a, b, m, plane[n].re[1][0], plane[n].im[1][0]);
		if (sp.harmonic)
		{
			write_map(dir + "E_nonlinear_" + suffix, comm[n], a, b, m, plane[n].re[0][1], plane[n].im[0][1]);
			if (dft_h)
				write_map(dir + "H_nonlinear_" + suffix, comm[n], a, b, m, plane[n].re[1][1], plane[n].im[1][1]);
		}
	};
	for (int m = 0; m < sp.nr_frequencies; ++m)
	{
		const std::string wl = std::to_string(std::lround((PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9))) + "nm";
		if (plane_on(2))
			full(2, m, wl);
		if (plane_on(1))
			full(1, m, wl);
		if (subdomain.plots_xz() && sp.graphene_optim)
			write_map(dir + "E_graphene_" + wl + "_XZ.h5", fdtd::grid::output_comm_plot_xz, 0, 2, m, graphene.re[0][0], graphene.im[0][0]);
		if (subdomain.plots_xz_2() && dft_xz_2)
			write_map(dir + "E_forward_" + wl + "_XZ.h5", fdtd::grid::output_comm_plot_xz_2, 0, 2, m, xz_2.re[0][0], xz_2.im[0][0]);
		if (subdomain.plots_xz_1() && dft_xz_1)
			write_map(dir + "E_backward_" + wl + "_XZ.h5", output_comm_plot_xz_1, 0, 2, m, xz_1.re[0][0], xz_1.im[0][0]);
		if (plane_on(0))
			full(0, m, wl);
	}
}
} // namespace fdtd::monitors
