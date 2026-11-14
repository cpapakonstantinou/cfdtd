/**
 * @file monitors/volume_dft.cc
 * @brief The volume DFT, monitors/volume_dft.h
 */
#include "monitors/volume_dft.h"
#include "monitors/state.h"
#include "core/views.h"
#include "grid/subdomain.h"
#include "fields/fields.h"
#include "materials/models.h"
#include "core/kernels.h"
#include <array>
#include <cmath>
#include "sources/state.h"
#include "monitors/spec.h"
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::fields::Hx;
using fdtd::fields::Hy;
using fdtd::fields::Hz;
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
/// Running DFT of a whole field over the inner region: one contiguous k-row kernel per (m, i, j)
void accumulate_volume_dft(double t, const std::array<view3<double>, 3> &f, const std::array<view4<double>, 3> &re, const std::array<view4<double>, 3> &im)
{
	const int nk = subdomain.inner_dim_z + 1;
	for (int m = 0; m < sp.nr_frequencies; ++m)
	{
		const double s = double(fdtd::sources::analysis_omega[m] * t);
		const double u = double(cos(s)), v = double(sin(s));
		for (int c = 0; c < 3; ++c)
			for (int i = 0; i <= subdomain.inner_dim_x; ++i)
				for (int j = 0; j <= subdomain.inner_dim_y; ++j)
					kernels::dft_accumulate(re[c][m].row(i, j), im[c][m].row(i, j), u, v,
											f[c].row(i + subdomain.inner_org_x, j + subdomain.inner_org_y) + subdomain.inner_org_z, nk);
	}
}
} // namespace

VolumeDft &volume_dft()
{
	static VolumeDft v;
	return v;
}

void VolumeDft::allocate()
{
	if (subdomain.cpml_only || !(sp.save_output || sp.abs_sigma || multipoles_wanted()))
		return;
	const int n = sp.nr_frequencies, nx = subdomain.inner_dim_x, ny = subdomain.inner_dim_y, nz = subdomain.inner_dim_z;
	auto real = [&](ctetrad<double> &a) {
		if (!a)
			a = field4<double>(n, nx, ny, nz);
	};
	for (int c = 0; c < 3; ++c)
	{
		for (int f = 0; f < 2; ++f)
			if (!field[f][c])
				field[f][c] = field4<std::complex<double>>(n, nx, ny, nz);
		real(re[0][c]);
		real(im[0][c]);
		if (multipoles_wanted() && (sp.shg || sp.thg) && n > 2)
		{
			real(re_pnl[c]);
			real(im_pnl[c]);
		}
		if (sp.save_output)
		{
			real(re[1][c]);
			real(im[1][c]);
		}
	}
}

void VolumeDft::release()
{
	for (int c = 0; c < 3; ++c)
	{
		for (int f = 0; f < 2; ++f)
		{
			re[f][c].reset();
			im[f][c].reset();
			field[f][c].reset();
		}
		re_pnl[c].reset();
		im_pnl[c].reset();
	}
}
void VolumeDft::sample_e(double t)
{
	accumulate_volume_dft(t, {Ex, Ey, Ez}, {re[0][0], re[0][1], re[0][2]}, {im[0][0], im[0][1], im[0][2]});
}

void VolumeDft::sample_h(double t)
{
	accumulate_volume_dft(t, {Hx, Hy, Hz}, {re[1][0], re[1][1], re[1][2]}, {im[1][0], im[1][1], im[1][2]});
}

void VolumeDft::normalize(int f)
{
	for (int c = 0; c < 3; ++c)
	{
		const view4<double> r = re[f][c], q = im[f][c];
		const view4<std::complex<double>> F = field[f][c];
		for (int m = 0; m < sp.nr_frequencies; ++m)
			for (int i = 0; i <= subdomain.inner_dim_x; ++i)
				for (int j = 0; j <= subdomain.inner_dim_y; ++j)
					for (int k = 0; k <= subdomain.inner_dim_z; ++k)
					{
						F[m][i][j][k] = std::complex<double>(r[m][i][j][k], q[m][i][j][k]) / fdtd::sources::dft_normalization[m];
						r[m][i][j][k] = F[m][i][j][k].real();
						q[m][i][j][k] = F[m][i][j][k].imag();
					}
	}
}

/// the same dft_normalization as the fields, the identity on the harmonic band, see PulseTransform, so
/// the contrast and the nonlinear parts of the total current add consistently
void VolumeDft::normalize_pnl()
{
	if (!(sp.shg || sp.thg) || sp.nr_frequencies <= 2)
		return;
	for (int c = 0; c < 3; ++c)
	{
		const view4<double> r = re_pnl[c], q = im_pnl[c];
		for (int m = sp.nr_frequencies / 2; m < sp.nr_frequencies; ++m)
			for (int i = 0; i < subdomain.inner_dim_x; ++i)
				for (int j = 0; j < subdomain.inner_dim_y; ++j)
					for (int k = 0; k < subdomain.inner_dim_z; ++k)
					{
						const std::complex<double> dft = std::complex<double>(r[m][i][j][k], q[m][i][j][k]) / fdtd::sources::dft_normalization[m];
						r[m][i][j][k] = dft.real();
						q[m][i][j][k] = dft.imag();
					}
	}
}

} // namespace fdtd::monitors
