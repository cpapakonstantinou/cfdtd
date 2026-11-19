/**
 * @file monitors/plane_dft.cc
 * @brief The running DFT on one plane and its normalisation, monitors/plane_dft.h
 *
 * Replaces updateDFT.cc, twelve entry points: E/H x field/harmonic x three orientations, and
 * DFTNormalizeXY / XZ / YZ.
 */
#include "monitors/plane_dft.h"
#include "monitors/state.h"
#include "fields/fields.h"
#include "core/views.h"
#include "grid/subdomain.h"
#include "grid/state.h"
#include "sources/state.h"
#include "core/kernels.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>
#include "monitors/spec.h"
using fdtd::grid::subdomain;

namespace
{
const fdtd::monitors::Spec &sp = fdtd::monitors::spec();        ///< what the monitors need from the run
} // namespace
using namespace fdtd::monitors; // transitional: module state, monitors/state.h

namespace fdtd::monitors
{
namespace
{
enum class Plane { YZ, XZ, XY };

/**
 * @brief Running DFT of up to three field components on one plane of the inner region
 * @verbatim
      for every frequency m and plane cell p:
          re[c][m][p] += f_c(p) cos(h w_m t)
          im[c][m][p] += f_c(p) sin(h w_m t)

      YZ plane, i fixed: p = (j, k)   accumulators [m][0][j][k]
      XZ plane, j fixed: p = (i, k)   accumulators [m][i][0][k]
      XY plane, k fixed: p = (i, j)   accumulators [m][i][j][0]
   @endverbatim
 * Each plane is gathered once into a contiguous buffer, then every frequency is one contiguous
 * kernels::dft_accumulate over the whole plane. Before, the field plane was read once per frequency,
 * and the XY plane with a stride of a whole k row per cell, 2-6x slower, see docs/dev/Roadmap.md.
 *
 * @param pad       plane position relative to the inner-region origin
 * @param harmonic  h: 1 for the fields, 2 for the nonlinear, second-harmonic, planes
 * @param f         field components x, y, z
 * @param re, im    accumulators of the components; a null `re` skips that component
 */
void plane_dft(double t, double harmonic, Plane plane, int pad, const std::array<view3<double>, 3> &f,
			  const std::array<view4<double>, 3> &re, const std::array<view4<double>, 3> &im)
{
	const int ni = subdomain.inner_dim_x + 1;
	const int nj = subdomain.inner_dim_y + 1;
	const int nk = subdomain.inner_dim_z + 1;
	const int i0 = subdomain.inner_org_x;
	const int j0 = subdomain.inner_org_y;
	const int k0 = subdomain.inner_org_z;
	// plane rows a and columns b, in the accumulators' order
	const int na = plane == Plane::YZ ? nj : ni;
	const int nb = plane == Plane::XY ? nj : nk;
	const std::size_t n = std::size_t(na) * nb;

	static std::vector<double> buffer;
	buffer.resize(3 * n);
	for (int c = 0; c < 3; ++c)
	{
		if (!re[c])
			continue;
		double *g = buffer.data() + c * n;
		for (int a = 0; a < na; ++a)
		{
			double *row = g + std::size_t(a) * nb;
			if (plane == Plane::YZ)
				std::copy_n(f[c].row(pad + i0, a + j0) + k0, nb, row);
			else if (plane == Plane::XZ)
				std::copy_n(f[c].row(a + i0, pad + j0) + k0, nb, row);
			else
				for (int b = 0; b < nb; ++b)
					row[b] = f[c][a + i0][b + j0][pad + k0];
		}
	}
	// the accumulator plane of each frequency must be one contiguous block of n values
	for (int c = 0; c < 3; ++c)
		if (re[c] && n > 1)
		{
			const view3<double> r = re[c][0];
			const double *last = plane == Plane::YZ ? r.row(0, na - 1) + nb - 1 : plane == Plane::XZ ? r.row(na - 1, 0) + nb - 1 : r.row(na - 1, nb - 1);
			if (last - r.row(0, 0) != std::ptrdiff_t(n) - 1)
				throw std::logic_error("plane_dft: accumulator plane is not contiguous");
		}

	for (int m = 0; m < sp.nr_frequencies; ++m)
	{
		const double s = double(harmonic * fdtd::sources::analysis_omega[m] * t);
		const double u = double(std::cos(s));
		const double v = double(std::sin(s));
		for (int c = 0; c < 3; ++c)
			if (re[c])
				kernels::dft_accumulate(re[c][m].row(0, 0), im[c][m].row(0, 0), u, v, buffer.data() + c * n, n);
	}
}

} // namespace

void sample_plane(int normal, int pad, int field, int set, double t, const array4<double> (&re)[3], const array4<double> (&im)[3])
{
	const fields::Fields &F = fields::fields();
	const int s = set == 1 && sp.harmonic ? 1 : 0;
	const array3<double> (&f)[3] = field == 0 ? F.E[s] : F.H[s];
	static const Plane kPlane[3] = {Plane::YZ, Plane::XZ, Plane::XY};
	plane_dft(t, set == 1 ? 2 : 1, kPlane[normal], pad, {f[0], f[1], f[2]}, {re[0], re[1], re[2]}, {im[0], im[1], im[2]});
}

void normalize_plane(int normal, const array4<double> (&re)[3], const array4<double> (&im)[3])
{
	const int ni = normal == 0 ? 0 : subdomain.inner_dim_x;
	const int nj = normal == 1 ? 0 : subdomain.inner_dim_y;
	const int nk = normal == 2 ? 0 : subdomain.inner_dim_z;
	for (int c = 0; c < 3; ++c)
	{
		const view4<double> r = re[c];
		const view4<double> q = im[c];
		if (!r || !q)
			continue;
		for (int m = 0; m < sp.nr_frequencies; ++m)
			for (int i = 0; i <= ni; ++i)
				for (int j = 0; j <= nj; ++j)
					for (int k = 0; k <= nk; ++k)
						if (phase_normalization)
						{
							const std::complex<double> dft = std::complex<double>(r[m][i][j][k], q[m][i][j][k]) / fdtd::sources::dft_normalization[m];
							r[m][i][j][k] = dft.real();
							q[m][i][j][k] = dft.imag();
						}
						else
						{
							r[m][i][j][k] /= fdtd::sources::amplitude_normalization[m];
							q[m][i][j][k] /= fdtd::sources::amplitude_normalization[m];
						}
	}
}
} // namespace fdtd::monitors
