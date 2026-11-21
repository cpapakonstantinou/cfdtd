/**
 * @file sources/waves.cc
 * @brief The plane-wave list of the superposition and Gaussian-beam sources, sources/waves.h
 *
 * Moved out of fdtd::app::initialise_run(), setupSimulation.cc, unchanged.
 */
#include "sources/waves.h"
#include "sources/state.h"
#include "io/text.h"
#include "sources/signal.h"
#include "materials/state.h"
#include "grid/state.h"
#include "grid/subdomain.h"
#include "grid/lattice.h"
#include "sources/spec.h"
using fdtd::constants::C0;
using fdtd::constants::PI;
using fdtd::constants::ROOTPROCESS;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::sources; // transitional: module state, sources/state.h

namespace fdtd::sources
{
bool build_waves(std::vector<PlaneWave> &waves, WaveFrame &frame)
{
	const auto &mp = fdtd::materials::model_parameters();
	waves.clear();
	if (!(sp.gauss_beam || sp.superposition))
		return true;

	// These are required for the plane wave superpositon formulas
	frame.focus[0] = (lat.cpml[0] + sp.centre[0]) * lat.d[0];
	frame.focus[1] = (lat.cpml[1] + sp.centre[1]) * lat.d[1];
	frame.focus[2] = (lat.cpml[2] + sp.centre[2]) * lat.d[2];

	const int i_tfsf_n = lat.cpml[0] + lat.tfsf_inset[0];
	const int j_tfsf_n = lat.cpml[1] + lat.tfsf_inset[1];
	const int k_tfsf_n = lat.cpml[2] + lat.tfsf_inset[2];

	const int nx_tfsf = lat.cells[0] - 2 * lat.cpml[0] - 2 * lat.tfsf_inset[0];
	const int ny_tfsf = lat.cells[1] - 2 * lat.cpml[1] - 2 * lat.tfsf_inset[1];
	const int nz_tfsf = lat.cells[2] - 2 * lat.cpml[2] - 2 * lat.tfsf_inset[2];

	const int i_tfsf_f = i_tfsf_n + nx_tfsf;
	const int j_tfsf_f = j_tfsf_n + ny_tfsf;
	const int k_tfsf_f = k_tfsf_n + nz_tfsf;

	frame.tfsf_n[0] = i_tfsf_n * lat.d[0];
	frame.tfsf_f[0] = i_tfsf_f * lat.d[0];

	frame.tfsf_n[1] = j_tfsf_n * lat.d[1];
	frame.tfsf_f[1] = j_tfsf_f * lat.d[1];

	frame.tfsf_n[2] = k_tfsf_n * lat.d[2];
	frame.tfsf_f[2] = k_tfsf_f * lat.d[2];

	// Construct vector of PlaneWave structs for the Gaussian beam
	if (sp.gauss_beam)
	{
		GaussianBeamData beamData;

		beamData.maxTheta = (PI / 180.) * sp.theta_max_degrees;
		beamData.nrTheta = sp.theta_intervals;
		beamData.deltaTheta = beamData.maxTheta / beamData.nrTheta;

		beamData.maxPhi = (PI / 180.) * sp.phi_max_degrees;
		beamData.nrPhi = sp.phi_intervals;
		beamData.deltaPhi = beamData.maxPhi / beamData.nrPhi;

		double fl_x = (sp.centre[0] - lat.tfsf_inset[0]) * lat.d[0];
		double fl_y = (sp.centre[1] - lat.tfsf_inset[1]) * lat.d[1];
		double fl_z = (sp.centre[2] - lat.tfsf_inset[2]) * lat.d[2];
		beamData.focalLength = sp.focal_factor * sqrt(fl_x * fl_x + fl_y * fl_y + fl_z * fl_z);
		beamData.beamWaist = sp.waist_factor * sin((PI / 180.) * sp.theta_max_degrees) * beamData.focalLength;

		if (fdtd::grid::world_rank == ROOTPROCESS)
		{
			std::string beamMessage = ("GAUSSIAN BEAM INFORMATION:\n\tWAIST:\t"
								  + std::to_string(beamData.beamWaist / 1.E-9)
								  + " nm\n\tFOCAL LENGTH:\t"
								  + std::to_string(beamData.focalLength / 1.E-9)
								  + " nm\n");
			fdtd::io::write_log(beamMessage, sp.reports_dir, sp.log_start, true);
		}

		if (sp.theta_intervals == 0 && sp.phi_intervals == 0)
		{
			bool useDerivative = false;
			double amp = 1.;
			// Set Offset / Delay for starttime of the wave
			frame.offset_delay = 0.0;
			build_oblique_plane_wave(beamData.maxTheta, beamData.maxPhi, amp, useDerivative, waves);
		}
		else
		{
			// Define delay for Gauss beam. Delay time is subtracted to every excitation at the TF/SF box
			frame.offset_delay = beamData.focalLength * mp.n_refr_index[sp.medium[0]] / C0;
			beamData.amplConst = ((mp.n_refr_index[sp.medium[0]] * beamData.focalLength / 2. / PI / C0)
								  * beamData.deltaTheta * beamData.deltaPhi);
			build_gaussian_plane_waves(beamData, waves);
		}
	}
	else if (sp.superposition)
	{
		// Set Offset / Delay for starttime of the wave
		frame.offset_delay = 0.0;
		// general plane wave superposition: the [wave.*] sections
		bool read_waves = fdtd::sources::build_plane_waves(plane_wave_input, waves, frame);
		if (!read_waves){
			return false;
		}
	}
	return true;
}
} // namespace fdtd::sources
