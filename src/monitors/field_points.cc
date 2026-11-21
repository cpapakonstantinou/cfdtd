/**
 * @file monitors/field_points.cc
 * @brief The point files, monitors/field_points.h
 *
 * Replaces updateFieldMonitors.cc, five copies of one writer, and its five flag globals.
 */
#include "monitors/field_points.h"
#include "io/text.h"
#include "monitors/state.h"
#include "core/views.h"
#include "grid/subdomain.h"
#include "materials/component.h"
#include "monitors/art_planes.h"
#include "monitors/dft_planes.h"
#include <cmath>
#include <cstdio>
#include <string>
#include "sources/state.h"
#include "monitors/spec.h"
#include "geometry/shapes.h"
using fdtd::constants::C0;
using fdtd::constants::PI;
using fdtd::constants::PI2;
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
/// one file: Ex, Ez of the E set of plane P at cell (i, j, k), for every frequency
void write_point(const std::string &name, const DftPlane &P, int i, int j, int k)
{
	FILE *f = fdtd::io::open_report(sp.reports_dir + "/" + name + ".txt", "lambda_nm\tabs_ex\tabs_ez\tabs_ex_over_ez\tphase_ex\tphase_ez\tphase_ex_minus_ez");
	if (!f)
		return;
	const view4<double> rex = P.re[0][0][0];
	const view4<double> iex = P.im[0][0][0];
	const view4<double> rez = P.re[0][0][2];
	const view4<double> iez = P.im[0][0][2];
	for (int m = 0; m < sp.nr_frequencies; ++m)
	{
		const double Ex_mod = std::sqrt(rex[m][i][j][k] * rex[m][i][j][k] + iex[m][i][j][k] * iex[m][i][j][k]);
		const double Ez_mod = std::sqrt(rez[m][i][j][k] * rez[m][i][j][k] + iez[m][i][j][k] * iez[m][i][j][k]);
		const double Ex_ph = std::atan2(iex[m][i][j][k], rex[m][i][j][k]) * 180 / PI;
		const double Ez_ph = std::atan2(iez[m][i][j][k], rez[m][i][j][k]) * 180 / PI;
		double Ex_Ez_ph = Ex_ph - Ez_ph;
		if (Ex_Ez_ph > 180)
			Ex_Ez_ph -= 360;
		if (Ex_Ez_ph < -180)
			Ex_Ez_ph += 360;
		std::fprintf(f, "%e\t %e\t %e\t %e\t %e\t %e\t %e\n", (PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9), Ex_mod, Ez_mod, Ex_mod / Ez_mod, Ex_ph,
					 Ez_ph, Ex_Ez_ph);
	}
	std::fclose(f);
}
} // namespace

void write_field_points()
{
	if (sp.grating)
		return;
	const bool antenna_family = (sp.shape == fdtd::geometry::shape::dipole) || (sp.shape == fdtd::geometry::shape::lateral_dipole) || (sp.shape == fdtd::geometry::shape::l_antenna) || (sp.shape == fdtd::geometry::shape::cylinder) || (sp.shape == fdtd::geometry::shape::x_shape) || (sp.shape == fdtd::geometry::shape::round_butterfly) ||
								(sp.shape == fdtd::geometry::shape::round_dipole) || (sp.shape == fdtd::geometry::shape::round_tripole) || (sp.shape == fdtd::geometry::shape::round_hybrid) || (sp.shape == fdtd::geometry::shape::one_side_round_dipole);
	const bool both_planes = subdomain.plots_xy() && subdomain.plots_yz();
	const int i = subdomain.i_animation;
	const int j = subdomain.j_animation;
	const int k = subdomain.k_animation;
	const ArtPlanes &art = art_planes();
	const DftPlanes &planes = dft_planes();

	if (both_planes && subdomain.art_far_y && sp.art && antenna_family)
		write_point("far_field_forward", art.plane[1], i, 0, k);
	if (both_planes && subdomain.art_nea_y && sp.art && antenna_family)
		write_point("far_field_backward", art.plane[0], i, 0, k);
	if (both_planes && subdomain.plots_xz_2() && (sp.shg || sp.thg) && (sp.shape == fdtd::geometry::shape::cylinder) && !sp.art)
		write_point("far_field_forward", planes.xz_2, i, 0, k);
	if (both_planes && subdomain.plots_xz_1() && (sp.shg || sp.thg) && (sp.shape == fdtd::geometry::shape::cylinder) && !sp.art)
		write_point("far_field_backward", planes.xz_1, i, 0, k);

	const bool near_shape = (sp.shape == fdtd::geometry::shape::dipole) || (sp.shape == fdtd::geometry::shape::lateral_dipole) || (sp.shape == fdtd::geometry::shape::x_shape) || (sp.shape == fdtd::geometry::shape::round_butterfly) || (sp.shape == fdtd::geometry::shape::round_dipole) ||
							(sp.shape == fdtd::geometry::shape::round_tripole) || (sp.shape == fdtd::geometry::shape::round_hybrid) || (sp.shape == fdtd::geometry::shape::one_side_round_dipole);
	if (both_planes && subdomain.plots_xz() && near_shape)
	{
		if (dft_plane_xz)
			write_point("near_field", planes.plane[1], i, 0, k);
		else if (dft_plane_yz)
			write_point("near_field", planes.plane[0], 0, j, k);
	}
}
} // namespace fdtd::monitors
