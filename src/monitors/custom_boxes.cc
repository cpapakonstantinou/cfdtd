/**
 * @file monitors/custom_boxes.cc
 * @brief The custom monitor boxes, monitors/custom_boxes.h
 *
 * Replaces readMonitors.cc, setupMonitors.cc, writeMonitors.cc, UpdateCustomMonitorDFTs and the
 * Monitor<T> struct, Monitors.h; its time instantiation carried unused DFT arrays.
 */
#include "monitors/custom_boxes.h"
#include "io/text.h"
#include "fields/fields.h"
#include "grid/subdomain.h"
#include "core/constants.h"
#include "materials/material_map.h"
#include "grid/state.h"
#include "materials/state.h"
#include "sources/state.h"
#include "io/binary.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include "monitors/spec.h"
using fdtd::constants::C0;
using fdtd::constants::PI2;
using fdtd::constants::ROOTPROCESS;
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
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
CustomBoxes &custom_boxes()
{
	static CustomBoxes c;
	return c;
}

namespace
{
CustomBox make_box(const MonitorBox &in)
{
	const int start[3] = {subdomain.shift_output_x, subdomain.shift_output_y, subdomain.shift_output_z};
	const int size[3] = {subdomain.inner_dim_x, subdomain.inner_dim_y, subdomain.inner_dim_z};
	const int center[3] = {in.center_x, in.center_y, in.center_z};
	const int sizes[3] = {in.size_x, in.size_y, in.size_z};

	CustomBox B;
	B.name = in.name;
	B.step = in.step;
	for (int a = 0; a < 3; a++)
	{
		B.global_min[a] = int(std::ceil(center[a] - 0.5 * sizes[a]));
		B.global_max[a] = int(std::ceil(center[a] + 0.5 * sizes[a]));
		B.global_count[a] = B.global_max[a] - B.global_min[a];
	}
	B.here = true;
	for (int a = 0; a < 3; a++)
		B.here = B.here && B.global_min[a] < start[a] + size[a] && B.global_max[a] >= start[a];
	if (B.here)
		for (int a = 0; a < 3; a++)
		{
			B.local_min[a] = B.global_min[a] < start[a] ? 0 : B.global_min[a] - start[a];
			B.local_max[a] = std::min(B.global_max[a] - start[a], size[a]);
			B.local_count[a] = B.local_max[a] - B.local_min[a];
			B.relative_start[a] = start[a] + B.local_min[a] - B.global_min[a];
		}

	// the communicator: ranks with inner cells and one of the box cells [global_min, global_max]
	bool member = true;
	for (int a = 0; a < 3; a++)
		member = member && size[a] > 0 && B.global_min[a] < start[a] + size[a] && B.global_max[a] >= start[a];
	int any = 0;
	int mine = member ? 1 : 0;
	MPI_Allreduce(&mine, &any, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
	if (!any && fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("monitors", "Custom Monitor outside of domain, ignoring...");
	MPI_Comm_split(MPI_COMM_WORLD, member ? 0 : MPI_UNDEFINED, fdtd::grid::world_rank, &B.comm);
	return B;
}
} // namespace

void CustomBoxes::release()
{
	for (std::vector<CustomBox> *set : {&frequency, &time})
	{
		for (CustomBox &B : *set)
			if (B.comm != MPI_COMM_NULL)
				MPI_Comm_free(&B.comm);
		set->clear();
	}
}

void CustomBoxes::build(const std::vector<MonitorBox> &boxes)
{
	release();
	for (const MonitorBox &in : boxes)
		(in.type == "frequency" ? frequency : time).push_back(make_box(in));
	for (CustomBox &B : frequency) // only where the rank holds a part, the old code also sized them, from unset counts, elsewhere
		if (B.here)
			for (int c = 0; c < 3; ++c)
				B.dft[c] = field4<std::complex<double>>(sp.nr_frequencies, B.local_count[0], B.local_count[1], B.local_count[2]);
}

void CustomBoxes::write_masks()
{
	auto write = [](CustomBox &B, const std::string &path) {
		// nrShift is local_min: the material maps have no CPML, the writer's PML flag just selects the shifted read
		save_as_binary_3d(path, matind_x, matind_y, matind_z, B.comm, B.relative_start, B.global_count, B.local_count, B.local_min,
					   false, true, true);
	};
	for (CustomBox &B : frequency)
		if (B.here)
			write(B, sp.fields_dir + "/E_DFT_" + B.name + "_Mask.h5");
	for (CustomBox &B : time)
		if (B.here)
			write(B, sp.frames_dir + "/E_Probe_" + B.name + "_Mask.h5");
}

void CustomBoxes::sample(double t, int step)
{
	const fields::Fields &F = fields::fields();
	for (CustomBox &B : frequency)
	{
		if (!B.here)
			continue;
		for (int m = 0; m < sp.nr_frequencies; m++)
		{
			const double s = double(fdtd::sources::analysis_omega[m] * t);
			const std::complex<double> phase = {double(std::cos(s)), double(std::sin(s))};
			for (int i = B.local_min[0]; i < B.local_max[0]; i++)
				for (int j = B.local_min[1]; j < B.local_max[1]; j++)
					for (int k = B.local_min[2]; k < B.local_max[2]; k++)
					{
						const int I = i + subdomain.inner_org_x;
						const int J = j + subdomain.inner_org_y;
						const int K = k + subdomain.inner_org_z;
						const int II = i - B.local_min[0];
						const int JJ = j - B.local_min[1];
						const int KK = k - B.local_min[2];
						for (int c = 0; c < 3; ++c)
							B.dft[c][m][II][JJ][KK] += F.E[0][c][I][J][K] * phase;
					}
		}
	}

	for (CustomBox &B : time)
	{
		if (B.here && step % B.step == 0)
		{
			MPI_Barrier(B.comm);
			const std::string path = sp.frames_dir + "/E_Probe_" + B.name + "_Frame_" + std::to_string(step / B.step) + ".h5";
			int shift[3] = {subdomain.i_cpml_min < 0 ? B.local_min[0] : subdomain.i_cpml_min + B.local_min[0],
							subdomain.j_cpml_min < 0 ? B.local_min[1] : subdomain.j_cpml_min + B.local_min[1],
							subdomain.k_cpml_min < 0 ? B.local_min[2] : subdomain.k_cpml_min + B.local_min[2]};
			save_as_binary_3d(path, Ex, Ey, Ez, B.comm, B.relative_start, B.global_count, B.local_count, shift, false, true, true);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

void CustomBoxes::write_dft()
{
	int shift[3] = {0, 0, 0};
	for (CustomBox &B : frequency)
		if (B.here)
			for (int m = 0; m < sp.nr_frequencies; m++)
			{
				const std::string wl = std::to_string(std::lround((PI2 * C0 / fdtd::sources::analysis_omega[m]) * double(1.e9)));
				save_as_binary_3d(sp.fields_dir + "/E_DFT_" + B.name + "_WL_" + wl + "nm.h5", B.dft[0][m], B.dft[1][m], B.dft[2][m],
							   B.comm, B.relative_start, B.global_count, B.local_count, shift, false, true, false);
			}
	MPI_Barrier(MPI_COMM_WORLD);
}
} // namespace fdtd::monitors
