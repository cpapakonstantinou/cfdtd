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
#include "monitors/exact_multipoles.h"
#include <new>
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
/// the 3D and plane field saves and the absorption integral read the whole volume
bool whole_volume() { return sp.save_output || sp.abs_sigma; }

/// the nonlinear polarisation of the full-current multipoles
bool pnl_wanted() { return multipoles_wanted() && (sp.shg || sp.thg) && sp.nr_frequencies > 2; }

/// Running DFT of a field over the box: one contiguous k-row kernel per (m, c, i, j)
void accumulate_volume_dft(double t, const std::array<view3<double>, 3> &f, array4<std::complex<double>> (&dft)[3], const CellBox &b)
{
	for (int m = 0; m < sp.nr_frequencies; ++m)
	{
		const double s = double(fdtd::sources::analysis_omega[m] * t);
		const double u = double(cos(s));
		const double v = double(sin(s));
		for (int c = 0; c < 3; ++c)
		{
			const view4<std::complex<double>> F = dft[c];
			for (int i = 0; i < b.n[0]; ++i)
				for (int j = 0; j < b.n[1]; ++j)
					kernels::dft_accumulate(F[m].row(i, j), u, v,
											f[c].row(i + b.lo[0] + subdomain.inner_org_x, j + b.lo[1] + subdomain.inner_org_y) + b.lo[2] + subdomain.inner_org_z, b.n[2]);
		}
	}
}

/// every value of a, frequency q + first, divided by S(w)
void divide_by_source(array4<std::complex<double>> &a, int first, int count, const CellBox &b)
{
	const view4<std::complex<double>> F = a;
	for (int q = 0; q < count; ++q)
		for (int i = 0; i < b.n[0]; ++i)
			for (int j = 0; j < b.n[1]; ++j)
				for (int k = 0; k < b.n[2]; ++k)
					F[q][i][j][k] /= fdtd::sources::dft_normalization[q + first];
}
} // namespace

VolumeDft &volume_dft()
{
	static VolumeDft v;
	return v;
}

bool VolumeDft::allocate()
{
	release();
	harmonic_start = sp.nr_frequencies / 2;
	if (subdomain.cpml_only || !(whole_volume() || multipoles_wanted()))
		return true;
	if (whole_volume())
		box = CellBox{{0, 0, 0}, {subdomain.inner_dim_x + 1, subdomain.inner_dim_y + 1, subdomain.inner_dim_z + 1}};
	else
		box = multipole_box();
	if (box.empty())
		return true;
	auto block = [&](int frequencies) {
		array4<std::complex<double>> a(frequencies, box.n[0], box.n[1], box.n[2]);
		a.zeros();
		return a;
	};
	try
	{
		for (int c = 0; c < 3; ++c)
		{
			field[0][c] = block(sp.nr_frequencies);
			if (sp.save_output)
				field[1][c] = block(sp.nr_frequencies);
			if (pnl_wanted())
				pnl[c] = block(sp.nr_frequencies - harmonic_start);
		}
	}
	catch (const std::bad_alloc &)
	{
		release();
		return false;
	}
	return true;
}

void VolumeDft::release()
{
	for (int c = 0; c < 3; ++c)
	{
		field[0][c].reset();
		field[1][c].reset();
		pnl[c].reset();
	}
	box = CellBox();
}

void VolumeDft::sample_e(double t)
{
	if (field[0][0])
		accumulate_volume_dft(t, {Ex, Ey, Ez}, field[0], box);
}

void VolumeDft::sample_h(double t)
{
	if (field[1][0])
		accumulate_volume_dft(t, {Hx, Hy, Hz}, field[1], box);
}

void VolumeDft::normalize(int f)
{
	for (int c = 0; c < 3; ++c)
		if (field[f][c])
			divide_by_source(field[f][c], 0, sp.nr_frequencies, box);
}

/// the same dft_normalization as the fields, the identity on the harmonic band, see PulseTransform, so
/// the contrast and the nonlinear parts of the total current add consistently
void VolumeDft::normalize_pnl()
{
	for (int c = 0; c < 3; ++c)
		if (pnl[c])
			divide_by_source(pnl[c], harmonic_start, sp.nr_frequencies - harmonic_start, box);
}

} // namespace fdtd::monitors
