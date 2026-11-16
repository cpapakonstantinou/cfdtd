/**
 * @file sources/tfsf_source.cc
 * @brief The TF/SF plane-wave source of this rank, sources/tfsf_source.h
 */
#include "sources/tfsf_source.h"
#include "sources/state.h"
#include "sources/signal.h"
#include "grid/subdomain.h"
#include "grid/state.h"
#include "sources/aux_grid.h"
#include <cmath>
#include "grid/lattice.h"
#include "sources/spec.h"
using fdtd::constants::INVVACIMP;
using fdtd::grid::subdomain;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::sources; // transitional: module state, sources/state.h

namespace fdtd::sources
{
TfsfSource TfsfSource::for_this_rank(const Signal &signal)
{
	TfsfSource pw;
	pw._signal = &signal;
	pw._active = subdomain.has_plane_wave && sp.tfsf && !sp.mode_source;
	if (!pw._active)
		return pw;
	const auto &S = subdomain;
	auto add = [](std::vector<Face> &list, bool on, Face f) {
		if (on)
			list.push_back(f);
	};
	add(pw._h_faces, S.pw_nea_y, Face::YNea);
	add(pw._e_faces, S.pw_nea_y, Face::YNea);
	add(pw._h_faces, S.pw_far_y, Face::YFar);
	add(pw._e_faces, S.pw_far_y, Face::YFar);
	if (!lat.periodic_x)
	{
		add(pw._h_faces, S.pw_nea_x && zpol, Face::XNea);
		add(pw._h_faces, S.pw_far_x && zpol, Face::XFar);
		add(pw._e_faces, S.pw_nea_x && xpol, Face::XNea);
		add(pw._e_faces, S.pw_far_x && xpol, Face::XFar);
	}
	if (!lat.periodic_z)
	{
		add(pw._h_faces, S.pw_nea_z && xpol, Face::ZNea);
		add(pw._h_faces, S.pw_far_z && xpol, Face::ZFar);
		add(pw._e_faces, S.pw_nea_z && zpol, Face::ZNea);
		add(pw._e_faces, S.pw_far_z && zpol, Face::ZFar);
	}
	return pw;
}

/**
 * Both polarisations at once, oblique, optionally circular: the z part of E (Hx) is shifted by circ_time_shift:
 *   Hx += sin(angle) / eta0 signal(t + shift),   Hz -= cos(angle) / eta0 signal(t)
 * One polarisation: Hx += signal(t) / eta0 and Hz -= signal(t) / eta0, the unused one stays zero in
 * the 3D grid: its faces are not in the lists.
 */
void TfsfSource::inject_h(const StepTime &time) const
{
	const double t = time.h;
	if (!_active)
		return;
	AuxGrid &g = aux_grid();
	const int lines = lat.periodic_x && lat.periodic_z ? 2 : 1;
	const int j = jpw_source_uni;
	if (xpol && zpol)
	{
		const double sx = sin(polarization_angle) * INVVACIMP;
		const double cz = cos(polarization_angle) * INVVACIMP;
		const double sig_x = _signal->value(t + sp.circ_time_shift);
		const double sig_z = _signal->value(t);
		for (int l = 0; l < lines; ++l)
		{
			g.line[l].h[1][j] += sx * sig_x;
			g.line[l].h[0][j] -= cz * sig_z;
		}
	}
	else
	{
		const double sig = INVVACIMP * _signal->value(t);
		for (int l = 0; l < lines; ++l)
		{
			g.line[l].h[1][j] += sig;
			g.line[l].h[0][j] -= sig;
		}
	}
	for (Face f : _h_faces)
		tfsf_h(f);
}

void TfsfSource::inject_e(const StepTime &) const
{
	for (Face f : _e_faces)
		tfsf_e(f);
}
} // namespace fdtd::sources
