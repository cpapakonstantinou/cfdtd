/**
 * @file sources/superposition.cc
 * @brief Plane-wave superposition source on the TF/SF faces, sources/superposition.h
 *
 * One face update for H and one for E serve the six faces, before: twelve functions in
 * gaussianBeam.cc; the face table and signs are in faces.h.
 */
#include "sources/superposition.h"
#include "sources/state.h"
#include "sources/signal.h"
#include "faces.h"
#include "materials/state.h"
#include "fields/fields.h"
#include "sources/spec.h"
using fdtd::constants::C0;
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
using namespace fdtd::sources; // transitional: module state, sources/state.h

namespace fdtd::sources
{
namespace
{
using namespace detail;

double direction_component(const PlaneWave &w, int c) { return c == 0 ? w.Dir_X : c == 1 ? w.Dir_Y : w.Dir_Z; }
double amp_e(const PlaneWave &w, int c) { return c == 0 ? w.AmpEx : c == 1 ? w.AmpEy : w.AmpEz; }
double amp_h(const PlaneWave &w, int c) { return c == 0 ? w.AmpHx : c == 1 ? w.AmpHy : w.AmpHz; }

/**
 * @brief Loops the cells of one face and calls cell(va, vb, To_a, Tm_a, To_b, Tm_b) with the projected
 *        positions of the integer (To) and half-integer in Tm points along a and b
 */
template <class Cell>
void for_face(const Box &X, const FaceAxes &F, const PlaneWave &w, Cell &&cell)
{
	for (int va = X.min[F.a]; va <= X.max[F.a]; ++va)
	{
		const int G = X.global(F.a, va);
		const double To_a = (G * X.d[F.a] - X.focus[F.a]) * direction_component(w, F.a), Tm_a = ((G + 0.5) * X.d[F.a] - X.focus[F.a]) * direction_component(w, F.a);
		for (int vb = X.min[F.b]; vb <= X.max[F.b]; ++vb)
		{
			const int K = X.global(F.b, vb);
			const double To_b = (K * X.d[F.b] - X.focus[F.b]) * direction_component(w, F.b), Tm_b = ((K + 0.5) * X.d[F.b] - X.focus[F.b]) * direction_component(w, F.b);
			cell(va, vb, To_a, Tm_a, To_b, Tm_b);
		}
	}
}

/// amplitude * signal, or its derivative, at t minus the travel time to the point with projections T[x, y, z]
double incident(const Signal &signal, const PlaneWave &w, double amp, double t, const double (&T)[3])
{
	const auto &mp = fdtd::materials::model_parameters();
	const double delay = (T[0] + T[1] + T[2]) * mp.n_refr_index[sp.medium[0]] / (C0 * w.dispersionFactor);
	return amp * (w.UseDerivative ? signal.derivative(t - delay) : signal.value(t - delay));
}

void face_h(const Signal &signal, const Box &X, const FaceAxes &F, bool far, const PlaneWave &w, double t)
{
	const view3<double> H[3] = {Hx, Hy, Hz};
	const double c = (far ? -F.s : F.s) * X.ch[F.n];
	const int p = far ? X.plane_f[F.n] : X.h_plane_n[F.n];
	const double T_n = ((far ? X.tfsf_f[F.n] : X.tfsf_n[F.n]) - X.focus[F.n]) * direction_component(w, F.n);
	for_face(X, F, w, [&](int va, int vb, double To_a, double Tm_a, double To_b, double Tm_b) {
		double Tb[3], Ta[3];
		Tb[F.n] = Ta[F.n] = T_n;
		Tb[F.a] = To_a, Tb[F.b] = Tm_b; // E_b at the half-integer point along b
		Ta[F.a] = Tm_a, Ta[F.b] = To_b; // E_a at the half-integer point along a
		const double e_b = incident(signal, w, amp_e(w, F.b), t, Tb);
		const double e_a = incident(signal, w, amp_e(w, F.a), t, Ta);
		if (vb != X.max[F.b])
			face_value(H[F.a], F, p, va, vb) += c * e_b;
		if (va != X.max[F.a])
			face_value(H[F.b], F, p, va, vb) -= c * e_a;
	});
}

void face_e(const Signal &signal, const Box &X, const FaceAxes &F, bool far, const PlaneWave &w, double t)
{
	const auto &mp = fdtd::materials::model_parameters();
	const view3<double> E[3] = {Ex, Ey, Ez};
	const double c = (far ? -F.s : F.s) * mp.ce_2[sp.medium[0]] * X.inv_d[F.n];
	const int p = far ? X.plane_f[F.n] : X.e_plane_n[F.n];
	// H half a cell outside the face
	const double T_n = ((far ? X.tfsf_f[F.n] + 0.5 * X.d[F.n] : X.tfsf_n[F.n] - 0.5 * X.d[F.n]) - X.focus[F.n]) * direction_component(w, F.n);
	for_face(X, F, w, [&](int va, int vb, double To_a, double Tm_a, double To_b, double Tm_b) {
		double Tb[3], Ta[3];
		Tb[F.n] = Ta[F.n] = T_n;
		Tb[F.a] = Tm_a, Tb[F.b] = To_b; // H_b at the half-integer point along a
		Ta[F.a] = To_a, Ta[F.b] = Tm_b; // H_a at the half-integer point along b
		const double h_b = incident(signal, w, amp_h(w, F.b), t, Tb);
		const double h_a = incident(signal, w, amp_h(w, F.a), t, Ta);
		if (va != X.max[F.a] && vb != 0)
			face_value(E[F.a], F, p, va, vb) -= c * h_b;
		if (vb != X.max[F.b] && va != 0)
			face_value(E[F.b], F, p, va, vb) += c * h_a;
	});
}

} // namespace

void superposition_h(const Signal &signal, const WaveFrame &frame, double t, const std::vector<PlaneWave> &waves)
{
	const Box X(frame);
	for (const PlaneWave &w : waves)
		for_rank_faces([&](const FaceAxes &F, bool far) { face_h(signal, X, F, far, w, t - frame.offset_delay); });
}

void superposition_e(const Signal &signal, const WaveFrame &frame, double t, const std::vector<PlaneWave> &waves)
{
	const Box X(frame);
	for (const PlaneWave &w : waves)
		for_rank_faces([&](const FaceAxes &F, bool far) { face_e(signal, X, F, far, w, t - frame.offset_delay); });
}
} // namespace fdtd::sources
