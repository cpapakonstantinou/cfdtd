/**
 * @file monitors/time_probes.cc
 * @brief The time-domain probes, monitors/time_probes.h
 *
 * Replaces SetupTimeDomainProbes, the PROBES table (legacy_parameters.h) and its global, and the
 * probe traces of the driver context.
 */
#include "monitors/time_probes.h"
#include "io/text.h"
#include "fields/fields.h"
#include "grid/subdomain.h"
#include "materials/material_map.h"
#include "grid/state.h"
#include "materials/state.h"
#include <cstdio>
#include <string>
#include "monitors/spec.h"
using fdtd::fields::index::y;
using fdtd::fields::index::z;
using fdtd::grid::subdomain;
using fdtd::materials::matind_x;
using fdtd::materials::matind_y;
using fdtd::materials::matind_z;

namespace
{
const fdtd::monitors::Spec &sp = fdtd::monitors::spec();        ///< what the monitors need from the run
} // namespace

namespace fdtd::monitors
{
TimeProbes &time_probes()
{
	static TimeProbes t;
	return t;
}

void TimeProbes::setup(const std::vector<ProbePoint> &probes)
{
	release();
	if (probes.empty() || subdomain.cpml_only)
		return;
	for (int i = subdomain.output_start[0]; i < subdomain.output_start[0] + subdomain.output_size[0] - 1; ++i)
		for (int j = subdomain.output_start[1]; j < subdomain.output_start[1] + subdomain.output_size[1] - 1; ++j)
			for (int k = subdomain.output_start[2]; k < subdomain.output_start[2] + subdomain.output_size[2] - 1; ++k)
				for (const ProbePoint &P : probes)
					if (i == P.x && j == P.y && k == P.z)
						point.push_back({int(i - subdomain.output_start[0]) + subdomain.inner_org_x, int(j - subdomain.output_start[1]) + subdomain.inner_org_y,
										 int(k - subdomain.output_start[2]) + subdomain.inner_org_z, P.x, P.y, P.z});
	if (point.empty())
		return;
	for (int c = 0; c < 3; ++c)
	{
		E[c] = array2<double>(int(point.size()), sp.steps + 1);
		E[c].zeros();
		H[c] = array2<double>(int(point.size()), sp.steps + 1);
		H[c].zeros();
	}
}

void TimeProbes::release()
{
	point.clear();
	for (int c = 0; c < 3; ++c)
	{
		E[c].reset();
		H[c].reset();
	}
}

void TimeProbes::sample(int step)
{
	const fields::Fields &F = fields::fields();
	for (size_t p = 0; p < point.size(); p++)
	{
		const ProbeCell &P = point[p];
		for (int c = 0; c < 3; ++c)
		{
			E[c][p][step] = F.E[0][c][P.x][P.y][P.z];
			H[c][p][step] = F.H[0][c][P.x][P.y][P.z];
		}
	}
}

void TimeProbes::write(int step) const
{
	for (size_t p = 0; p < point.size(); p++)
	{
		const ProbeCell &P = point[p];
		const std::string path = sp.frames_dir + "/probe_" + std::to_string(P.gx) + "_" + std::to_string(P.gy) + "_" + std::to_string(P.gz) + ".txt";
		FILE *f = fdtd::io::open_report(path, "");
		if (!f)
			continue;
		const int x = P.x - subdomain.inner_org_x;
		const int y = P.y - subdomain.inner_org_y;
		const int z = P.z - subdomain.inner_org_z;
		std::fprintf(f, "# material slots: %d %d %d\n", matind_x[x][y][z], matind_y[x][y][z], matind_z[x][y][z]);
		std::fprintf(f, "# ex\tey\tez\thx\thy\thz\n");
		for (int s = 0; s < step; s++)
			std::fprintf(f, "%e\t%e\t%e\t%e\t%e\t%e\n", E[0](p, s), E[1](p, s), E[2](p, s), H[0](p, s), H[1](p, s), H[2](p, s));
		std::fclose(f);
	}
}
} // namespace fdtd::monitors
