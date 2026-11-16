/**
 * @file monitors/art_planes.cc
 * @brief The reflection and transmission planes of periodic runs, monitors/art_planes.h
 */
#include "monitors/art_planes.h"
#include "monitors/state.h"
#include "monitors/plane_dft.h"
#include "grid/subdomain.h"
#include "core/views.h"
#include "grid/state.h"
#include "grid/lattice.h"
#include "monitors/spec.h"
using fdtd::constants::I_UNIT;
using fdtd::grid::subdomain;

namespace
{
const fdtd::monitors::Spec &sp = fdtd::monitors::spec();        ///< what the monitors need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::monitors; // transitional: module state, monitors/state.h

namespace fdtd::monitors
{
ArtPlanes &art_planes()
{
	static ArtPlanes a;
	return a;
}

bool ArtPlanes::has(int side) const
{
	return sp.art && !subdomain.cpml_only && (side == 0 ? subdomain.art_nea_y : subdomain.art_far_y);
}

void ArtPlanes::allocate()
{
	for (int side = 0; side < 2; ++side)
		if (sp.art && (side == 0 ? subdomain.art_nea_y : subdomain.art_far_y))
			for (int f = 0; f < 2; ++f)
				for (int c : {0, 2})
				{
					plane[side].re[f][0][c] = field4<double>(sp.nr_frequencies, subdomain.inner_dim_x, 0, subdomain.inner_dim_z);
					plane[side].im[f][0][c] = field4<double>(sp.nr_frequencies, subdomain.inner_dim_x, 0, subdomain.inner_dim_z);
				}
}

void ArtPlanes::release()
{
	for (DftPlane &P : plane)
		for (int f = 0; f < 2; ++f)
			for (int c = 0; c < 3; ++c)
			{
				P.re[f][0][c].reset();
				P.im[f][0][c].reset();
			}
}

/// as before: the near plane E then H, the far plane H then E
void ArtPlanes::sample(double t_e, double t_h)
{
	auto e = [&](DftPlane &P, int j) { sample_plane(1, j, 0, 0, t_e, P.re[0][0], P.im[0][0]); };
	auto h = [&](DftPlane &P, int j) { sample_plane(1, j, 1, 0, t_h, P.re[1][0], P.im[1][0]); };
	if (has(0))
	{
		e(plane[0], subdomain.j_art_1);
		h(plane[0], subdomain.j_art_1);
	}
	if (has(1))
	{
		h(plane[1], subdomain.j_art_2);
		e(plane[1], subdomain.j_art_2);
	}
}

void ArtPlanes::normalize()
{
	for (int side = 0; side < 2; ++side)
		if (has(side))
			for (int f = 0; f < 2; ++f)
				normalize_plane(1, plane[side].re[f][0], plane[side].im[f][0]);
}

double ArtPlanes::flux(int side, double out[]) const
{
	const DftPlane &P = plane[side];
	const view4<double> reEx = P.re[0][0][0];
	const view4<double> imEx = P.im[0][0][0];
	const view4<double> reEz = P.re[0][0][2];
	const view4<double> imEz = P.im[0][0][2];
	const view4<double> reHx = P.re[1][0][0];
	const view4<double> imHx = P.im[1][0][0];
	const view4<double> reHz = P.re[1][0][2];
	const view4<double> imHz = P.im[1][0][2];
	double q_normalisation_area = 0;
	for (int m = 0; m < sp.nr_frequencies; ++m)
	{
		q_normalisation_area = double(0.);
		double P_f = double(0.);
		for (int i = 0; i < subdomain.inner_dim_x; ++i)
		{
			const double dsurf1 = lat.d[0];
			for (int k = 0; k < subdomain.inner_dim_z; ++k)
			{
				const double dsurf2 = lat.d[2];
				q_normalisation_area += dsurf1 * dsurf2;
				const std::complex<double> cEx = reEx[m][i][0][k] + I_UNIT * imEx[m][i][0][k];
				const std::complex<double> cEz = reEz[m][i][0][k] + I_UNIT * imEz[m][i][0][k];
				const std::complex<double> ccHx = reHx[m][i][0][k] - I_UNIT * imHx[m][i][0][k];
				const std::complex<double> ccHz = reHz[m][i][0][k] - I_UNIT * imHz[m][i][0][k];
				const std::complex<double> zz = cEz * ccHx - cEx * ccHz;
				P_f += zz.real() * dsurf1 * dsurf2;
			}
		}
		out[m] = P_f;
	}
	return q_normalisation_area;
}
} // namespace fdtd::monitors
