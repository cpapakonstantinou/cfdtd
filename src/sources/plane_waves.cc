/**
 * @file sources/plane_waves.cc
 * @author FDTD-Lab, https://github.com/FDTD-Lab
 * @version 0.1
 * @date 2022-09-23
 *
 * @brief Plane waves of a wave_superposition source, [wave.*] sections
 *
 * @copyright Copyright 2022
 */
#include "sources/types.h"
#include "io/text.h"
#include "sources/state.h"
#include "sources/signal.h"
#include "materials/state.h"
#include "grid/state.h"
#include "grid/subdomain.h"
#include "io/binary.h"
#include <math.h>
#include <string>
#include <vector>
#include "grid/lattice.h"
#include "sources/spec.h"
#include "sources/kinds.h"
using fdtd::constants::C0;
using fdtd::constants::PI2;
using fdtd::constants::ROOTPROCESS;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace
using namespace fdtd::sources; // transitional: module state, sources/state.h

// Define Helper function
std::vector<double> get_corner_first_incidence(const WaveFrame &, double, double, double);
// function calculating the numerical dispersion factor using wavevector
double get_dispersion_factor(double, double, double);
double root_dispersion_relation(double, double, double, double);
/// function defining the derivative of the numerical dispersion relation
double derivative_root_dispersion_relation(double, double, double, double);

bool fdtd::sources::build_plane_waves(const std::vector<PlaneWave> &input, std::vector<PlaneWave> &returnVector, WaveFrame &frame)
{
	const auto &mp = fdtd::materials::model_parameters();
	const bool correctDispersion = plane_wave_dispersion;
	for (size_t i = 0; i < input.size(); i++)
	{
		returnVector.push_back(input[i]);
		// Correct the magnetic field amplitudes by refractive index of the first medium
		returnVector[i].AmpHx *= mp.n_refr_index[sp.medium[0]];
		returnVector[i].AmpHy *= mp.n_refr_index[sp.medium[0]];
		returnVector[i].AmpHz *= mp.n_refr_index[sp.medium[0]];

		// identify corner of first incidence
		std::vector<double> r0 = get_corner_first_incidence(frame, returnVector[i].Dir_X,
														 returnVector[i].Dir_Y,
														 returnVector[i].Dir_Z);

		// calculate the temporary offset delay
		double tmpDelay = (r0[0] * returnVector[i].Dir_X) +
						  (r0[1] * returnVector[i].Dir_Y) +
						  (r0[2] * returnVector[i].Dir_Z);
		tmpDelay = sqrt(tmpDelay * tmpDelay) * mp.n_refr_index[sp.medium[0]] / C0;

		// check if offset delay is larger than current and set if true
		frame.offset_delay = (tmpDelay > frame.offset_delay) ? tmpDelay : frame.offset_delay;

		// Calculate the numerical dispersion factor for this plane wave
		if (correctDispersion)
			returnVector[i].dispersionFactor = get_dispersion_factor(
				returnVector[i].Dir_X, returnVector[i].Dir_Y,
				returnVector[i].Dir_Z);
		else
			returnVector[i].dispersionFactor = 1.0;
	}

	return true;
}

/// @brief Function to return a vector from the center of the TFSF box to the
/// corner where a plane wave is first incident
/// @param kx x-component of direction vector
/// @param ky y-component of direction vector
/// @param kz z-component of direction vector
/// @return three dimensional vector pointing from center to the corner
std::vector<double> get_corner_first_incidence(const WaveFrame &frame, double kx, double ky, double kz)
{
	// initialize vector
	std::vector<double> corner;
	// check for x component
	double tmp = (kx >= 0) ? frame.tfsf_n[0] - frame.focus[0] : frame.tfsf_f[0] - frame.focus[0];
	corner.push_back(tmp);
	// y component
	tmp = (ky >= 0) ? frame.tfsf_n[1] - frame.focus[1] : frame.tfsf_f[1] - frame.focus[1];
	corner.push_back(tmp);
	// z component
	tmp = (kz >= 0) ? frame.tfsf_n[2] - frame.focus[2] : frame.tfsf_f[2] - frame.focus[2];
	corner.push_back(tmp);
	// return vector
	return corner;
}

/// @brief Calculates the dispersion factor associated with a plane wave
/// @details This implements the Newton Raphson method to calculate the
/// numerical wave vector. From this is it obtains the dipsersion factor, i.e.
/// the factor by which the speed of light gets smaller in the numerical grid.
/// @param kx
/// @param ky
/// @param kz
/// @return dispersion factor
double get_dispersion_factor(double kx, double ky, double kz)
{
	// if the source is not Gaussian we can not correct for dispersion
	if (sp.signal_type != pulse::gaussian)
		return 1.0;

	// Define the current wave vector
	double k = PI2 * run_signal().carrier_frequency() / C0;
	double analyticalRoot = root_dispersion_relation(kx, ky, kz, k);
	// current best guess for numerical wavevector
	double kNumerical = k;
	// counter for the number of iterations
	int counter = 0;
	const int maxIterations = 5000;
	// Implement Newton Raphson to find the numerical wave vector
	while (std::abs(root_dispersion_relation(kx, ky, kz, kNumerical)) > 1e-5 &&
		   counter < maxIterations)
	{
		// update numerical wave vector
		kNumerical =
			kNumerical -
			root_dispersion_relation(kx, ky, kz, kNumerical) /
				derivative_root_dispersion_relation(kx, ky, kz, kNumerical);
		// increase counter
		counter++;
	}
	// define the root value that we have found.
	double currRoot = root_dispersion_relation(kx, ky, kz, kNumerical);

	// check if counter exceeded
	if (counter == maxIterations)
	{
		// is analytical wave vector better than numerically calculated
		if (std::abs(analyticalRoot) < std::abs(currRoot))
		{
			if (fdtd::grid::world_rank == ROOTPROCESS)
				fdtd::io::log("sources", "Could not find numerical wave vector. Using Analytical "
					   "expression! \n");
			// Return 1.0 as dispersion factor, i.e. no correction for
			// dispersion
			return 1.0;
		}
	}
	// return dispersion factor
	return PI2 * run_signal().carrier_frequency() / (C0 * kNumerical);
}

/// @brief Root function for the dispersion relation. Returns zero if numerical
/// dispersion relation is matched.
/// @details This function defines the numerical dispersion relation for a given
/// direction of the plane wave and a given value for the absolute of the
/// numerical wavevector. If the function returns zero the numerical dispersion
/// relation is satisfied.
/// @param kx x-component of unit direction vector
/// @param ky y-component of unit direction vector
/// @param kz z-component of unit direction vector
/// @param r Absolute of the numerical wave vector
/// @return difference between LHS and RHS of numerical dispersion relation
double root_dispersion_relation(double kx, double ky, double kz, double r)
{
	// Use a temporary variable for better readability
	double tmp = -pow(1 / (C0 * sp.dt) * sin(PI2 * run_signal().carrier_frequency() * sp.dt / 2), 2.0);
	tmp += pow(1 / lat.d[0] * sin(r * kx * lat.d[0] / 2), 2.0);
	tmp += pow(1 / lat.d[1] * sin(r * ky * lat.d[1] / 2), 2.0);
	tmp += pow(1 / lat.d[2] * sin(r * kz * lat.d[2] / 2), 2.0);
	// Return value of tmp as it defines the dispersion relation
	return tmp;
}

/// @brief Derivative of the root_dispersion_relation() function
/// @param kx x-component of the unit direction vector
/// @param ky y-component of the unit direction vector
/// @param kz z-component of the unit direction vector
/// @param r Absolute of the numerical wave vector
/// @return derivative of root_dispersion_relation()
double derivative_root_dispersion_relation(double kx, double ky, double kz,
										double r)
{
	// Use temporary variable for calculating the derivative.
	double tmp = kx / lat.d[0] * sin(r * kx * lat.d[0] / 2.) * cos(r * kx * lat.d[0] / 2.);
	tmp += ky / lat.d[1] * sin(r * ky * lat.d[1] / 2.) * cos(r * ky * lat.d[1] / 2.);
	tmp += kz / lat.d[2] * sin(r * kz * lat.d[2] / 2.) * cos(r * kz * lat.d[2] / 2.);
	return tmp;
}