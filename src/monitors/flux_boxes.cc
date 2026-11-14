/**
 * @file monitors/flux_boxes.cc
 * @brief The DFT planes of the two flux boxes, monitors/flux_boxes.h
 *
 * One loop over boxes, faces and field sets replaces the 192 allocations, 192 releases and 96
 * sampling calls written out per array before.
 */
#include "monitors/flux_boxes.h"
#include "monitors/state.h"
#include "grid/subdomain.h"
#include "grid/state.h"
#include "monitors/plane_dft.h"
#include "monitors/spec.h"
using fdtd::grid::subdomain;

namespace
{
const fdtd::monitors::Spec &sp = fdtd::monitors::spec();        ///< what the monitors need from the run
} // namespace
using namespace fdtd::monitors; // transitional: module state, monitors/state.h

namespace fdtd::monitors
{
FluxBoxes &flux_boxes()
{
	static FluxBoxes f;
	return f;
}

namespace
{
/// does this rank hold face (n, side) of box b
bool has_box_face(int b, int n, int side)
{
	const auto &S = subdomain;
	const bool f[2][3][2] = {{{S.qsca_nea_x, S.qsca_far_x}, {S.qsca_nea_y, S.qsca_far_y}, {S.qsca_nea_z, S.qsca_far_z}},
							 {{S.qtot_nea_x, S.qtot_far_x}, {S.qtot_nea_y, S.qtot_far_y}, {S.qtot_nea_z, S.qtot_far_z}}};
	return f[b][n][side];
}

/// local plane index of face (n, side) of box b
int plane_index(int b, int n, int side)
{
	const int p[2][3][2] = {{{subdomain.i_qsca_nea, subdomain.i_qsca_far}, {subdomain.j_qsca_nea, subdomain.j_qsca_far}, {subdomain.k_qsca_nea, subdomain.k_qsca_far}},
							{{subdomain.i_qtot_nea, subdomain.i_qtot_far}, {subdomain.j_qtot_nea, subdomain.j_qtot_far}, {subdomain.k_qtot_nea, subdomain.k_qtot_far}}};
	return p[b][n][side];
}

} // namespace

void FluxBoxes::allocate()
{
	if (!q_coeff)
		return;
	const int sets = sp.hydrodynamics ? 2 : 1;
	for (int b = 0; b < 2; ++b)
		for (int n = 0; n < 3; ++n)
		{
			if (!has_box_face(b, n, 0) && !has_box_face(b, n, 1))
				continue;
			const int dims[3] = {n == 0 ? 0 : subdomain.inner_dim_x, n == 1 ? 0 : subdomain.inner_dim_y, n == 2 ? 0 : subdomain.inner_dim_z};
			for (int side = 0; side < 2; ++side) // both faces, as before, when the rank has either
				for (int f = 0; f < 2; ++f)
					for (int s = 0; s < sets; ++s)
						for (int c = 0; c < 3; ++c)
							if (c != n)
							{
								box[b].face[n][side].re[f][s][c] = field4<double>(sp.nr_frequencies, dims[0], dims[1], dims[2]);
								box[b].face[n][side].im[f][s][c] = field4<double>(sp.nr_frequencies, dims[0], dims[1], dims[2]);
							}
		}
}

void FluxBoxes::release()
{
	for (FluxBox &B : box)
		for (auto &axis : B.face)
			for (FluxFace &F : axis)
				for (int f = 0; f < 2; ++f)
					for (int s = 0; s < 2; ++s)
						for (int c = 0; c < 3; ++c)
						{
							F.re[f][s][c].reset();
							F.im[f][s][c].reset();
						}
}

/**
 * Order as before: box scattering then total; faces x near, x far, y near, y far, z near, z far;
 * per face H, E, then, hydrodynamics, harmonic E, harmonic H.
 */
void FluxBoxes::sample(double t_e, double t_h)
{
	auto run = [&](FluxFace &F, int f, int s, int n, int p) {
		const double tt = f == 0 ? t_e : t_h;
		sample_plane(n, p, f, s, tt, F.re[f][s], F.im[f][s]);
	};
	for (int b = 0; b < 2; ++b)
		for (int n = 0; n < 3; ++n)
			for (int side = 0; side < 2; ++side)
			{
				if (!has_box_face(b, n, side))
					continue;
				FluxFace &F = box[b].face[n][side];
				const int p = plane_index(b, n, side);
				run(F, 1, 0, n, p);
				run(F, 0, 0, n, p);
				if (sp.hydrodynamics)
				{
					run(F, 0, 1, n, p);
					run(F, 1, 1, n, p);
				}
			}
}

void FluxBoxes::normalize()
{
	for (int b = 0; b < 2; ++b)
		for (int n = 0; n < 3; ++n)
			if (has_box_face(b, n, 0) || has_box_face(b, n, 1))
				for (int side = 0; side < 2; ++side)
					for (int f = 0; f < 2; ++f)
						normalize_plane(n, box[b].face[n][side].re[f][0], box[b].face[n][side].im[f][0]);
}
} // namespace fdtd::monitors
