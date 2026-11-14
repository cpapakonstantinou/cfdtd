/**
 * @file sources/faces.h
 * @brief Internal to the sources module: the six faces of the TF/SF box for the analytic sources
 *
 * A face has a normal n and tangential axes a < b. The correction signs on the near face:
 * @verbatim
      H_a += s ch_n E_b_inc     H_b -= s ch_n E_a_inc
      E_a -= s c_n  H_b_inc     E_b += s c_n  H_a_inc      (c_n = ce2 / d_n)
   @endverbatim
 * with s = +1 for the y face (a = x, b = z) and -1 for the x (y, z) and z (x, y) faces; the far face
 * has -s. Faces run y, x, z near, then y, x, z far, edge cells get two corrections in that order.
 */
#ifndef FDTD_SOURCES_FACES_H
#define FDTD_SOURCES_FACES_H

#include "grid/subdomain.h"
#include "core/views.h"
#include "grid/state.h"
#include "update/state.h"
#include "sources/state.h"
#include "sources/spec.h"
#include "grid/lattice.h"

namespace fdtd::sources::detail
{
struct FaceAxes
{
	int n, a, b; ///< normal, tangential axes (a < b)
	double s;    ///< sign of the H_a correction on the near face
};
inline constexpr FaceAxes kFaces[3] = {{1, 0, 2, 1}, {0, 1, 2, -1}, {2, 0, 1, -1}}; // y, x, z

/// the TF/SF box of this rank per axis
struct Box
{
	int min[3], max[3];  ///< local range, inclusive
	long org[3];         ///< global index = local - org
	double d[3], inv_d[3], tfsf_n[3], tfsf_f[3], ch[3];
	double focus[3]; ///< focus of the superposition, Gaussian beam
	int h_plane_n[3], e_plane_n[3], plane_f[3]; ///< the face planes, local index along n

	explicit Box(const WaveFrame &frame)
	{
		using fdtd::grid::subdomain;
		using namespace fdtd::materials; // ipw_near_h, jpw_near_h, kpw_near_h
		const auto &S = subdomain;
		const int mn[3] = {subdomain.ipw_min, subdomain.jpw_min, subdomain.kpw_min}, mx[3] = {subdomain.ipw_max, subdomain.jpw_max, subdomain.kpw_max};
		const long og[3] = {S.inner_org_x - fdtd::grid::lattice().cpml[0] - subdomain.output_start[0], S.inner_org_y - fdtd::grid::lattice().cpml[1] - subdomain.output_start[1],
							S.inner_org_z - fdtd::grid::lattice().cpml[2] - subdomain.output_start[2]};
		const double dd[3] = {fdtd::grid::lattice().d[0], fdtd::grid::lattice().d[1], fdtd::grid::lattice().d[2]}, id[3] = {fdtd::grid::lattice().inv_d[0], fdtd::grid::lattice().inv_d[1], fdtd::grid::lattice().inv_d[2]}, tn[3] = {frame.tfsf_n[0], frame.tfsf_n[1], frame.tfsf_n[2]},
					 tf[3] = {frame.tfsf_f[0], frame.tfsf_f[1], frame.tfsf_f[2]}, c[3] = {fdtd::sources::chx, fdtd::sources::chy, fdtd::sources::chz}, fo[3] = {frame.focus[0], frame.focus[1], frame.focus[2]};
		const int hp[3] = {ipw_near_h, jpw_near_h, kpw_near_h}, ep[3] = {subdomain.ipw_1, subdomain.jpw_1, subdomain.kpw_1}, fp[3] = {subdomain.ipw_2, subdomain.jpw_2, subdomain.kpw_2};
		for (int q = 0; q < 3; ++q)
		{
			min[q] = mn[q], max[q] = mx[q], org[q] = og[q], d[q] = dd[q], inv_d[q] = id[q], tfsf_n[q] = tn[q], tfsf_f[q] = tf[q];
			ch[q] = c[q], focus[q] = fo[q], h_plane_n[q] = hp[q], e_plane_n[q] = ep[q], plane_f[q] = fp[q];
		}
	}
	int global(int q, int local) const { return int(local - org[q]); }
};

/// cell, n at plane p, a at va, b at vb, of f
inline double &face_value(view3<double> f, const FaceAxes &F, int p, int va, int vb)
{
	int x[3];
	x[F.n] = p;
	x[F.a] = va;
	x[F.b] = vb;
	return f[x[0]][x[1]][x[2]];
}

/// face(F, far) for this rank's faces, y, x, z near then y, x, z far
template <class Face>
void for_rank_faces(Face &&face)
{
	const auto &S = fdtd::grid::subdomain;
	const bool nea[3] = {S.pw_nea_y, S.pw_nea_x, S.pw_nea_z}, far[3] = {S.pw_far_y, S.pw_far_x, S.pw_far_z};
	for (int f = 0; f < 3; ++f)
		if (nea[f])
			face(kFaces[f], false);
	for (int f = 0; f < 3; ++f)
		if (far[f])
			face(kFaces[f], true);
}
} // namespace fdtd::sources::detail

#endif
