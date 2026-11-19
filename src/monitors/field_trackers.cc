/**
 * @file monitors/field_trackers.cc
 * @brief The field trackers, monitors/field_trackers.h
 *
 * Replaces GetMultipolePositions, setupSources.cc, the tracker block of the step DFT and of the
 * checkpoint output, and 18 globals, the point lists, counts and the twelve per-component arrays.
 */
#include "monitors/field_trackers.h"
#include "io/text.h"
#include "monitors/state.h"
#include "fields/fields.h"
#include "grid/subdomain.h"
#include "grid/state.h"
#include "sources/state.h"
#include <cmath>
#include <cstdio>
#include <string>
#include "monitors/spec.h"
using fdtd::constants::C0;
using fdtd::constants::PI;
using fdtd::constants::PI2;
using fdtd::constants::ROOTPROCESS;
using fdtd::grid::subdomain;

namespace
{
const fdtd::monitors::Spec &sp = fdtd::monitors::spec();        ///< what the monitors need from the run
} // namespace
using namespace fdtd::monitors; // transitional: module state, monitors/state.h

namespace fdtd::monitors
{
FieldTrackers &field_trackers()
{
	static FieldTrackers f;
	return f;
}

void FieldTrackers::release()
{
	point.clear();
	for (auto &set : dft)
		for (auto &field : set)
			for (auto &component : field)
				component.reset();
}

void FieldTrackers::setup(int center_x, int center_y, int center_z)
{
	release();
	const int R = sphere_rad_multipole;
	const int Ntheta = ntheta_multipole;
	const int Nphi = nphi_multipole;
	std::vector<int> X(Ntheta * Nphi);
	std::vector<int> Y(Ntheta * Nphi);
	std::vector<int> Z(Ntheta * Nphi);
	for (int i = 0, n = 0; i < Ntheta; i++)
		for (int j = 0; j < Nphi; j++, n++)
		{
			const double phi = PI2 / Nphi * j;
			const double theta = PI / Ntheta * i;
			X[n] = int(R * std::sin(theta) * std::cos(phi));
			Y[n] = int(R * std::sin(theta) * std::sin(phi));
			Z[n] = int(R * std::cos(theta));
		}

	point.clear();
	for (int i = subdomain.output_start[0]; i < subdomain.output_start[0] + subdomain.output_size[0] - 1; ++i)
		for (int j = subdomain.output_start[1]; j < subdomain.output_start[1] + subdomain.output_size[1] - 1; ++j)
			for (int k = subdomain.output_start[2]; k < subdomain.output_start[2] + subdomain.output_size[2] - 1; ++k)
				for (int n = 0; n < Ntheta * Nphi; n++)
					if (i - center_x == X[n] && j - center_y == Y[n] && k - center_z == Z[n])
						point.push_back({int(i - subdomain.output_start[0]) + subdomain.inner_org_x, int(j - subdomain.output_start[1]) + subdomain.inner_org_y,
										 int(k - subdomain.output_start[2]) + subdomain.inner_org_z, n});

	int here = int(point.size());
	int total = 0;
	MPI_Reduce(&here, &total, 1, MPI_INT, MPI_SUM, ROOTPROCESS, MPI_COMM_WORLD);
	if (fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("monitors", "TOTAL POINTS = %d\n", total);

	if (point.empty())
		return;
	const int sets = sp.harmonic ? 2 : 1;
	for (int s = 0; s < sets; ++s)
		for (int f = 0; f < 2; ++f)
			for (int c = 0; c < 3; ++c)
				dft[s][f][c] = field3<double>(int(point.size()), sp.nr_frequencies, 2);
}

void FieldTrackers::sample(double t_e, double t_h)
{
	const fields::Fields &F = fields::fields();
	const int sets = sp.harmonic ? 2 : 1;
	for (size_t p = 0; p < point.size(); p++)
	{
		const TrackPoint &P = point[p];
		for (int m = 0; m < sp.nr_frequencies; ++m)
			for (int s = 0; s < sets; ++s)
			{
				const double w = s == 0 ? fdtd::sources::analysis_omega[m] : 2.0 * fdtd::sources::analysis_omega[m];
				for (int c = 0; c < 3; ++c)
				{
					const double e = F.E[s][c][P.x][P.y][P.z];
					const double h = F.H[s][c][P.x][P.y][P.z];
					dft[s][0][c][p][m][0] += std::cos(w * t_e) * e;
					dft[s][1][c][p][m][0] += std::cos(w * t_h) * h;
					dft[s][0][c][p][m][1] += std::sin(w * t_e) * e;
					dft[s][1][c][p][m][1] += std::sin(w * t_h) * h;
				}
			}
	}
}

void FieldTrackers::write() const
{
	static const char *name[2] = {"Efield_Tracker_", "Efield_NL_Tracker_"};
	const int sets = sp.harmonic ? 2 : 1;
	for (size_t p = 0; p < point.size(); p++)
		for (int s = 0; s < sets; ++s)
		{
			const TrackPoint &P = point[p];
			const std::string path = sp.frames_dir + "/" + name[s] + std::to_string(P.index) + ".txt";
			FILE *f = std::fopen(path.c_str(), "w");
			if (!f)
				continue;
			std::fprintf(f, "%ld\t%ld\t%ld\n", P.x + subdomain.output_start[0] - subdomain.inner_org_x - sp.centre[0],
						 P.y + subdomain.output_start[1] - subdomain.inner_org_y - sp.centre[1], P.z + subdomain.output_start[2] - subdomain.inner_org_z - sp.centre[2]);
			const double harmonic = s == 0 ? 1.0 : 2.0;
			const array3<double>(&a)[2][3] = dft[s];
			for (int m = 0; m < sp.nr_frequencies; ++m)
				std::fprintf(f, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", C0 / (harmonic * fdtd::sources::analysis_omega[m] / (2 * PI)),
							 a[0][0][p][m][0], a[0][0][p][m][1], a[0][1][p][m][0], a[0][1][p][m][1], a[0][2][p][m][0], a[0][2][p][m][1],
							 a[1][0][p][m][0], a[1][0][p][m][1], a[1][1][p][m][0], a[1][1][p][m][1], a[1][2][p][m][0], a[1][2][p][m][1]);
			std::fclose(f);
		}
}
} // namespace fdtd::monitors
