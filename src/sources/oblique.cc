/**
 * @file sources/oblique.cc
 * @brief Plane waves through a layered background on the TF/SF faces, sources/oblique.h
 *
 * One face update for H and one for E serve the six faces, before: twelve functions in
 * obliquePlaneWave.cc; the face table and signs are in faces.h. The incident field of component c
 * is fields[c] (E) or fields[3 + c] (H) of OneDimensionalSimulator::return_fields at its Yee point:
 * @verbatim
      H face corrections:  E_a at (a + 1/2, b, n),  E_b at (a, b + 1/2, n), n on the face
      E face corrections:  H_b at (a + 1/2, b, n -+ 1/2),  H_a at (a, b + 1/2, n -+ 1/2), half a cell outside
   @endverbatim
 * E coefficient: ce_2 of the first, near y face, or second, far y face, medium, or of the layer at
 * the cell's y on the x and z faces. The 1D simulator is queried in the order of the previous code.
 */
#include "sources/oblique.h"
#include "sources/one_dim_simulator.h"
#include "faces.h"
#include "materials/state.h"
#include "fields/fields.h"
#include "sources/spec.h"
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::fields::Hx;
using fdtd::fields::Hy;
using fdtd::fields::Hz;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
} // namespace

namespace fdtd::sources
{
namespace
{
using namespace detail;

/// incident fields at the point with coordinates p[x, y, z]
std::array<double, 6> incident_fields(OneDimensionalSimulator &ods, const double (&p)[3], double t) { return ods.return_fields(p[0], p[1], p[2], t); }

/// calls cell(va, vb, x) for the face cells, x the integer-point coordinates, x[n] = face
template <class Cell>
void for_face(const Box &X, const FaceAxes &F, double face, Cell &&cell)
{
	for (int va = X.min[F.a]; va <= X.max[F.a]; ++va)
		for (int vb = X.min[F.b]; vb <= X.max[F.b]; ++vb)
		{
			double x[3];
			x[F.n] = face;
			x[F.a] = X.d[F.a] * X.global(F.a, va);
			x[F.b] = X.d[F.b] * X.global(F.b, vb);
			cell(va, vb, x);
		}
}
} // namespace

void oblique_h(OneDimensionalSimulator &ods, double t)
{
	const Box X(ods.frame());
	const view3<double> H[3] = {Hx, Hy, Hz};
	for_rank_faces([&](const FaceAxes &F, bool far) {
		const double c = (far ? -F.s : F.s) * X.ch[F.n];
		const int p = far ? X.plane_f[F.n] : X.h_plane_n[F.n];
		for_face(X, F, far ? X.tfsf_f[F.n] : X.tfsf_n[F.n], [&](int va, int vb, const double (&x)[3]) {
			double pa[3] = {x[0], x[1], x[2]};
			double pb[3] = {x[0], x[1], x[2]};
			pa[F.a] = x[F.a] + 0.5 * X.d[F.a]; // E_a
			pb[F.b] = x[F.b] + 0.5 * X.d[F.b]; // E_b
			const double e_a = incident_fields(ods, pa, t)[F.a];
			if (va != X.max[F.a])
				face_value(H[F.b], F, p, va, vb) -= c * e_a;
			const double e_b = incident_fields(ods, pb, t)[F.b];
			if (vb != X.max[F.b])
				face_value(H[F.a], F, p, va, vb) += c * e_b;
		});
	});
}

void oblique_e(OneDimensionalSimulator &ods, double t)
{
	const auto &mp = fdtd::materials::model_parameters();
	const Box X(ods.frame());
	const view3<double> E[3] = {Ex, Ey, Ez};
	for_rank_faces([&](const FaceAxes &F, bool far) {
		const double s = far ? -F.s : F.s;
		const int p = far ? X.plane_f[F.n] : X.e_plane_n[F.n];
		const double face = far ? X.tfsf_f[F.n] : X.tfsf_n[F.n];
		for_face(X, F, face, [&](int va, int vb, const double (&x)[3]) {
			double pb[3] = {x[0], x[1], x[2]};
			double pa[3] = {x[0], x[1], x[2]};
			pb[F.n] = pa[F.n] = far ? face + 0.5 * X.d[F.n] : face - 0.5 * X.d[F.n];
			pb[F.a] = x[F.a] + 0.5 * X.d[F.a]; // H_b
			pa[F.b] = x[F.b] + 0.5 * X.d[F.b]; // H_a
			const std::array<double, 6> h_b = incident_fields(ods, pb, t);
			const int m = F.n == 1 ? (far ? sp.medium[2] : sp.medium[0]) : ods.get_medium_index(x[1]);
			const double c = s * mp.ce_2[m] * X.inv_d[F.n];
			if (va != X.max[F.a] && vb != 0)
				face_value(E[F.a], F, p, va, vb) -= c * h_b[3 + F.b];
			const double h_a = incident_fields(ods, pa, t)[3 + F.a];
			if (vb != X.max[F.b] && va != 0)
				face_value(E[F.b], F, p, va, vb) += c * h_a;
		});
	});
}
} // namespace fdtd::sources
