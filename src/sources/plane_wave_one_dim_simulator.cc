#include "sources/types.h"
#include "io/text.h"
#include "sources/signal.h"
#include "core/constants.h"
#include "materials/state.h"
#include "sources/plane_wave_one_dim_simulator.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>
#include "grid/lattice.h"
#include "sources/spec.h"
using fdtd::constants::C0;
using fdtd::constants::EP0;
using fdtd::constants::MU0;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice
} // namespace

bool PlaneWaveOneDimensionalSimulator::init(PlaneWave &tmpWave, const WaveFrame &frame)
{
	_frame = frame;
	_signal = &fdtd::sources::run_signal();
	// set plane wave struct
	pw = tmpWave;

	set_geometry_parameters(pw);

	set_grid_parameters();

	set_velocity_and_timestep();

	set_medium_indices();

	initialize_field_vectors();

	// run the simulation
	return run_simulation();
}

bool PlaneWaveOneDimensionalSimulator::set_geometry_parameters(PlaneWave &tmpWave)
{
	if (tmpWave.Dir_Y < 0.0)
		fdtd::io::log("sources",  "WARNING: You use a plane wave that is propagating in negative y-direction. "
				"The implementation assumes a wave propagating in postive y-direction. Hence"
				"this might lead to problems! \n"); 

	// calculate mx and mz
	if (std::abs(tmpWave.Dir_X) < 1e-5 && std::abs(tmpWave.Dir_Z) < 1e-5)
	{
		// default to m-vector in x-direction if PW only in y direction
		mx = 1.0;
		mz = 0.0;
	} else {
		
		double denominator =
			sqrt(tmpWave.Dir_X * tmpWave.Dir_X + tmpWave.Dir_Z * tmpWave.Dir_Z);
		mx = tmpWave.Dir_X / denominator;
		mz = tmpWave.Dir_Z / denominator;
	}

	// incident angle with y axis
	theta = acos(tmpWave.Dir_Y);

	// set the initial update fields
	initialCartesianFields = {tmpWave.AmpEx, tmpWave.AmpEy, tmpWave.AmpEz,
								tmpWave.AmpHx, tmpWave.AmpHy, tmpWave.AmpHz};
	initialUpdateFields =
		cartesian_fields_to_update_fields(initialCartesianFields);

	// set flag for derivative signal
	useDerivativeSignal = tmpWave.UseDerivative;

	return true;
}

/**
 * @brief sets the parameters read from the input file for this 1D simulation 
*/
bool PlaneWaveOneDimensionalSimulator::set_grid_parameters()
{
	gridLength = 3 * lat.cells[1];
	gridstep = lat.d[1];
	gridShiftIdx = lat.cells[1];
	numCpmlCells = 3 * lat.cpml[1];
	cpmlNearEndIdx = numCpmlCells;
	cpmlFarStartIdx = gridLength - numCpmlCells;
	tfsfStartIdx = lat.tfsf_inset[1] + lat.cpml[1] + gridShiftIdx - 1;

	return true;
}

bool PlaneWaveOneDimensionalSimulator::set_velocity_and_timestep()
{
	const auto &mp = fdtd::materials::model_parameters();
	// tangential phase velocity 
	tangentialPhaseVelocity = C0 / (mp.n_refr_index[sp.medium[0]] * sin(theta));

	// Determining the maximum propagation velocity in the 3 media
	maxPropagationVelocity = 1 / sqrt(EP0 * MU0 * (mp.epre[sp.medium[0]] - mp.epre[sp.medium[0]] * sin(theta) * sin(theta)));
	maxPropagationVelocity = std::max(maxPropagationVelocity, 1 / sqrt(EP0 * MU0 * (mp.epre[sp.medium[1]] - mp.epre[sp.medium[0]] * sin(theta) * sin(theta))));
	maxPropagationVelocity = std::max(maxPropagationVelocity, 1 / sqrt(EP0 * MU0 * (mp.epre[sp.medium[2]] - mp.epre[sp.medium[0]] * sin(theta) * sin(theta))));

	// Determining a time step to satisfy the Courant criterion
	timestep = sp.dt; 
	double k = 1.0;
	courantNumber = maxPropagationVelocity * timestep / gridstep;
	while(courantNumber >= 0.5){
		timestep = sp.dt / (2.0 * k + 1.0);
		courantNumber = maxPropagationVelocity * timestep / gridstep;
		k = k + 1;
	}

	// Time Step counters
	numInitialTimesteps = sp.one_dim_steps * (k-1);
	currTimestep = 0;

	// time delay for this individual plane wave to synchronise all plane waves at the domain center
	timeDelay = (mx * _frame.focus[0] + mz * _frame.focus[2]) / tangentialPhaseVelocity;

	std::array<double, 3> distCornerCenter =  get_distance_center_to_corner_first_incidence();
	timeDelay -= (distCornerCenter[1] * pw.Dir_Y) * mp.n_refr_index[sp.medium[0]] / C0;

	return true;
}

bool PlaneWaveOneDimensionalSimulator::set_medium_indices()
{
	const auto &mp = fdtd::materials::model_parameters();
	// define where the substrate starts <-- taken from CreateBackground functionality
	int deepY = sp.centre[1] - ceil(sp.middle_layer / lat.d[1]);
	// set the relative permittivities for the 1D simulation
	for ( int j = 0; j < gridLength; j++)
	{
		// create a shifted j, that matches to the grid of 3D domain
		int shiftedJ = j - gridShiftIdx - lat.cpml[1];
		// this is similar to the material initialization in CreateBackground
		if (shiftedJ < deepY)
		{
			relativePermittivity.push_back(mp.epre[sp.medium[0]]);
			mediumIndices.push_back(sp.medium[0]);
		}
		else if (shiftedJ >= deepY && shiftedJ < sp.centre[1])
		{
			relativePermittivity.push_back(mp.epre[sp.medium[1]]);
			mediumIndices.push_back(sp.medium[1]);
		}	
		else if (shiftedJ >= sp.centre[1])
		{
			relativePermittivity.push_back(mp.epre[sp.medium[2]]);
			mediumIndices.push_back(sp.medium[2]);
		}	
	}
	return true;
}

bool PlaneWaveOneDimensionalSimulator::initialize_field_vectors()
{
	// Initialize empty vectors for zeroth time step
	std::vector<double> initZeros(gridLength, 0.0);
	Eh.push_back(initZeros);
	Hh.push_back(initZeros);
	Hy.push_back(initZeros);
	Ee.push_back(initZeros);
	He.push_back(initZeros);
	Ey.push_back(initZeros);
	times.push_back(0.0);

	// resize the auxiliary fields for CPML regions
	psiEhFar.resize(numCpmlCells);
	psiHhFar.resize(numCpmlCells);
	psiEeFar.resize(numCpmlCells);
	psiHeFar.resize(numCpmlCells);
	psiEhNear.resize(numCpmlCells);
	psiHhNear.resize(numCpmlCells);
	psiEeNear.resize(numCpmlCells);
	psiHeNear.resize(numCpmlCells);
	
	// fill the auxiliary field with zeros as initial values
	std::fill(psiEhFar.begin(), psiEhFar.end(), 0.0);
	std::fill(psiHhFar.begin(), psiHhFar.end(), 0.0);
	std::fill(psiEeFar.begin(), psiEeFar.end(), 0.0);
	std::fill(psiHeFar.begin(), psiHeFar.end(), 0.0);
	std::fill(psiEhNear.begin(), psiEhNear.end(), 0.0);
	std::fill(psiHhNear.begin(), psiHhNear.end(), 0.0);
	std::fill(psiEeNear.begin(), psiEeNear.end(), 0.0);
	std::fill(psiHeNear.begin(), psiHeNear.end(), 0.0);

	return true;
}

int PlaneWaveOneDimensionalSimulator::get_medium_index(double y)
{
	// Set correct index in alignment to the 3D grid
	int j = std::floor(y / gridstep) + gridShiftIdx; 
	return mediumIndices[j];
}

double PlaneWaveOneDimensionalSimulator::signal_function(double t)
{
	if( useDerivativeSignal)
	{
		return _signal->derivative(t - _frame.offset_delay);
	} else
	{
		return _signal->value(t - _frame.offset_delay);
	}
}

/**
 * @brief Transforms cartesian fields, electric/magnetic, into the coordinate
 * system used for the 1D simulation
 * @param fields given as {Ex, Ey, Ez, Hx, Hy, Hz}
 * @return fields in new coordinate system {E^h, H^h, hy, E^e, H^e, ey}
 */ 
std::array<double, 6>
PlaneWaveOneDimensionalSimulator::cartesian_fields_to_update_fields(
	const std::array<double, 6> &fields)
{
	// Electric fields
	double em = mx * fields[0] + mz * fields[2];
	double en = mz * fields[0] - mx * fields[2];
	double ey = fields[1];
	// Magnetic fields
	double hm = mx * fields[3] + mz * fields[5];
	double hn = mz * fields[3] - mx * fields[5];
	double hy = fields[4];
	// return the fields
	return {en, -1.0 * hm, hy, em, hn, ey};
}

/**
 * @brief Transforms the fields from the coordinate system used in the 1D
 * simulation to cartesian coordinates
 * @param fields holds the fields used for the 1D simulation {E^h, H^h, H_y,E^e, H^e, ey}
 * @return fields in cartesian coordinate system {ex, ey, ez, hx, H_y, hz}
 */
std::array<double, 6>
PlaneWaveOneDimensionalSimulator::update_fields_to_cartesian_fields(
	const std::array<double, 6> &fields)
{
	// electric fields
	double ex = (mx * fields[3] + mz * fields[0]);
	double ez = (mz * fields[3] - mx * fields[0]);
	double ey = fields[5];
	// magnetic fields
	double hx = (mx * -1.0 * fields[1] + mz * fields[4]);
	double hz = (mz * -1.0 * fields[1] - mx * fields[4]);
	double hy = fields[2];
	// build new fields array
	return {ex, ey, ez, hx, hy, hz};
}

/**
 * @brief Updates fields to the next time step
 * @return Boolean if update proceeded correctly
*/ 
bool PlaneWaveOneDimensionalSimulator::propagate_timestep()
{
	// check that length are all the same
	int n = Eh.size();
	assert(Hh.size() == size_t(n));
	assert(Hy.size() == size_t(n));

	// set references to previous time step
	std::vector<double> hhPrev = Hh.back();
	std::vector<double> ehPrev = Eh.back();
	std::vector<double> hyPrev = Hy.back();
	std::vector<double> hePrev = He.back();
	std::vector<double> eePrev = Ee.back();
	std::vector<double> eyPrev = Ey.back();

	auto [hhNew, heNew] = update_magnetic_fields_in_domain(hhPrev, hePrev, ehPrev, eePrev);

	update_magnetic_fields_in_cpml_region(hhNew, heNew, ehPrev, hhPrev, eePrev, hePrev);

	auto [ehNew, eeNew] = update_electric_fields_in_domain(hhNew, heNew, ehPrev, eePrev);

	update_electric_fields_in_cpml_region(ehNew, eeNew, hhNew, ehPrev, heNew, eePrev);

	auto [hyNew, eyNew] = update_fields_ycomponents_in_domain(ehPrev, heNew, hyPrev, eyPrev);

	// Attach the new vectors to the storage of all vectors
	Hh.push_back(hhNew);
	Eh.push_back(ehNew);
	Hy.push_back(hyNew);
	He.push_back(heNew);
	Ee.push_back(eeNew);
	Ey.push_back(eyNew);

	// update current timestep
	currTimestep += 1;
	times.push_back(currTimestep * timestep);

	return true;
}

/**
 * @brief Runs the simulator that is setup in the init function
 * @return bool if simulation finished.
*/ 
bool PlaneWaveOneDimensionalSimulator::run_simulation()
{
	bool check = true;
	for (int n = 0; n < numInitialTimesteps; n++)
	{
		check = check && propagate_timestep();
	}

	return check;
}

/**
 * @brief Calculates the fields in cartesian fields at the given coordinates
 * @param x location in x-direction
 * @param y location in y-direction
 * @param z location in z-direction
 * @param t time at which to get the fields
 * @return fields at the given location and time
*/
std::array<double, 6> PlaneWaveOneDimensionalSimulator::return_fields(double x, double y, double z, double t)
{
	// calculate coordinates in update field system
	double m = mx * x + mz * z;

	// define time at which to evaluate the fields
	double tCorr = t - m / tangentialPhaseVelocity + timeDelay;

	// check that required time is already simulated
	while(tCorr > times.back() - 3 * timestep){
		// propagate time step
		propagate_timestep();

		// remove the front of the vectors
		Eh.erase(Eh.begin());
		Hh.erase(Hh.begin());
		Hy.erase(Hy.begin());
		Ee.erase(Ee.begin());
		He.erase(He.begin());
		Ey.erase(Ey.begin());
		times.erase(times.begin());
	}

	// check the sizes of the array
	assert(Eh.size() == size_t(numInitialTimesteps + 1));
	assert(Hh.size() == size_t(numInitialTimesteps + 1));
	assert(Hy.size() == size_t(numInitialTimesteps + 1));
	assert(Ee.size() == size_t(numInitialTimesteps + 1));
	assert(He.size() == size_t(numInitialTimesteps + 1));
	assert(Ey.size() == size_t(numInitialTimesteps + 1));
	assert(times.size() == size_t(numInitialTimesteps + 1));

	// indices surrounding the time step
	int n = floor( tCorr / timestep) - floor( times.front() / timestep); 

	// return zero if time is before simulated
	if (n <= 0)
		return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	// TODO: This should not be necessary
	if (times[n] > tCorr)
		n = n - 1;

	// Check that indices are correct
	assert(times[n] <= tCorr + 1e-10);
	assert(times[n+1] >= tCorr - 1e-10);

	// time index for magnetic fields
	int nMag = floor(tCorr / timestep + 0.5) - floor(times.front() / timestep + 0.5);

	// Check that time index is correct
	assert(times[nMag] - 0.5 * timestep <= tCorr + 1e-10);
	assert(times[nMag] + 0.5 * timestep >= tCorr - 1e-10);

	// Spatial index for fields at integers (Eh, Ee, Hy)
	int j = floor(y / gridstep);
	assert(j * gridstep <= y);
	assert((j + 1) * gridstep >= y);

	// spatial index for half integer fields (Hh, He, Ey)
	int jStag = floor(y / gridstep + 0.5);
	assert((jStag - 0.5) * gridstep <= y + 1e-10);
	assert((jStag + 0.5) * gridstep >= y - 1e-10);
	
	// factors for linear time interpolation of electric field
	double ct1 = (times[n + 1] - tCorr) / (times[n + 1] - times[n]);
	double ct2 = (tCorr - times[n]) / (times[n + 1] - times[n]);

	// factor for linear time interpolation of magnetic fields
	double ct1Mag = (times[nMag] + 0.5 * timestep - tCorr) / timestep;
	double ct2Mag = (tCorr - (times[nMag] - 0.5 * timestep)) / timestep;

	// factors for linear spatial interpolation of integer fields
	double cy1 = ((j+1) * gridstep - y) / gridstep;
	double cy2 = (y - j * gridstep) / gridstep;

	// factors for linear spatial interpolation of half-integer fields
	double cy1Stag = ((jStag + 0.5) * gridstep - y) / gridstep;
	double cy2Stag = (y - (jStag - 0.5) * gridstep) / gridstep;

	// shift the spatial indices according to 
	j = j + gridShiftIdx;
	jStag = jStag + gridShiftIdx;

	// interpolations of Eh and Ee
	double ehTmp = ct1 * (cy1 * Eh[n][j] + cy2 * Eh[n][j + 1]) + ct2 * (cy1 * Eh[n+1][j] + cy2 * Eh[n+1][j + 1]);
	double eeTmp = ct1 * (cy1 * Ee[n][j] + cy2 * Ee[n][j + 1]) + ct2 * (cy1 * Ee[n+1][j] + cy2 * Ee[n+1][j + 1]);
	
	// interpolation of Ey
	double eyTmp = ct1 * (cy1Stag * Ey[n][jStag] + cy2Stag * Ey[n][jStag + 1]) + ct2 * (cy1Stag * Ey[n+1][jStag] + cy2Stag * Ey[n+1][jStag + 1]);

	// interpolation of Hh and He
	double hhTmp = ct1Mag * (cy1Stag * Hh[nMag][jStag] + cy2Stag * Hh[nMag][jStag + 1]) + ct2Mag * (cy1Stag * Hh[nMag + 1][jStag] + cy2Stag * Hh[nMag + 1][jStag + 1]);
	double heTmp = ct1Mag * (cy1Stag * He[nMag][jStag] + cy2Stag * He[nMag][jStag + 1]) + ct2Mag * (cy1Stag * He[nMag + 1][jStag] + cy2Stag * He[nMag + 1][jStag + 1]);

	// interpolation of Hy
	double hyTmp = ct1Mag * (cy1 * Hy[nMag][j] + cy2 * Hy[nMag][j+1]) + ct2Mag * (cy1 * Hy[nMag + 1][j] + cy2 * Hy[nMag + 1][j + 1]);

	// convert fields to cartesian
	std::array<double, 6> fields = update_fields_to_cartesian_fields({ehTmp, hhTmp, hyTmp, eeTmp, heTmp, eyTmp});

	// return cartesian fields 
	return fields;
}

/**
 * @brief Returns the coefficients required for field updates in CPML regions
 * @param Spatial index in y direction at which the cpml coefficients are required
 * @return Coefficients b, c, cb and k
*/
std::array<double, 3> PlaneWaveOneDimensionalSimulator::get_cpml_coefficients_at_index(int j)
{
	double k, sigmaMax, sigma, alpha, b, c;

	// Depth at which the current update occurs in the CPML region
	double cpmlPenetrationDepth;
	if (j <= cpmlNearEndIdx)
		cpmlPenetrationDepth = (double) cpmlNearEndIdx - (double) j;
	else 
		cpmlPenetrationDepth = (double) j - (double) cpmlFarStartIdx;

	double mScaling = 3.5;
	double kMax = 10.0;
	double alphaMax = 0.3;
	double ma = 1.0;

	k = 1.0 + (kMax - 1.0) * pow( cpmlPenetrationDepth / numCpmlCells, mScaling);
	sigmaMax = 0.15 * 0.8 * (mScaling + 1.0) / (gridstep * sqrt(MU0 / EP0 * relativePermittivity[j]));
	sigma = sigmaMax * pow( cpmlPenetrationDepth / numCpmlCells, mScaling);
	alpha = alphaMax * pow((numCpmlCells - cpmlPenetrationDepth) / numCpmlCells, ma);
	b = exp(-1.0 * (sigma / EP0 / k + alpha / EP0) * timestep);
	c = sigma / (k * (sigma + alpha * k)) * (b - 1.0);

	// return only the relevant variables
	return {b, c, k};
}

/**
 *  @brief Updates the Magnetic Fields inside the domain region, not CPMLs
 * @param hhPrev - vector containing the Hh fields at the previous timestep
 * @param hePrev - vector containing the He fields at the previous timestep
 * @param ehPrev - vector containing the Eh fields at the previous timestep
 * @param eePrev - vector containing the Ee fields at the previous timestep
 * @return Vectors of fields Hh and He at next timestep
*/
std::array<std::vector<double>, 2> PlaneWaveOneDimensionalSimulator::update_magnetic_fields_in_domain(std::vector<double>& hhPrev, std::vector<double>& hePrev, std::vector<double>& ehPrev, std::vector<double>& eePrev)
{
	// setup vectors for new fields
	std::vector<double> hhNew(gridLength, 0.0);
	std::vector<double> heNew(gridLength, 0.0);

	// Loop through domain region
	for (int i = cpmlNearEndIdx; i < cpmlFarStartIdx + 1 ; i++)
	{
		// Update equation for Hh
		hhNew[i] = hhPrev[i] - timestep / (MU0 * gridstep) * (ehPrev[i] - ehPrev[i-1]);

		// Update equation for He
		double prefactorHe = timestep * relativePermittivity[i] / (MU0 * gridstep * (relativePermittivity[i] - relativePermittivity[0] * sin(theta) * sin(theta)));
		heNew[i] = hePrev[i] - prefactorHe * (eePrev[i] - eePrev[i-1]);

		// TFSF correction terms
		if (i == tfsfStartIdx)
		{
			double signalTime = timestep * currTimestep;
			hhNew[i] += timestep / (MU0 * gridstep) * initialUpdateFields[0] * signal_function(signalTime);
			heNew[i] += prefactorHe * initialUpdateFields[3] * signal_function(signalTime);
		}
	}
	return {hhNew, heNew};
}

/**
 *  @brief Updates the Electric Fields inside the domain region, not CPMLs
 * @param hhNew - vector containing the Hh fields at the next timestep
 * @param heNew - vector containing the He fields at the next timestep
 * @param ehPrev - vector containing the Eh fields at the previous timestep
 * @param eePrev - vector containing the Ee fields at the previous timestep
 * @return Vectors of fields Eh and Ee at next timestep
*/
std::array<std::vector<double>, 2> PlaneWaveOneDimensionalSimulator::update_electric_fields_in_domain(std::vector<double>& hhNew, std::vector<double>& heNew, std::vector<double>& ehPrev, std::vector<double>& eePrev)
{
	// setup vectors for fields at next time step
	std::vector<double> ehNew(gridLength, 0.0);
	std::vector<double> eeNew(gridLength, 0.0);

	// Loop through the domain region
	for ( int i = cpmlNearEndIdx; i < cpmlFarStartIdx + 1; i++)
	{
		// Update equation for Eh
		double prefactorEh = timestep / (gridstep * EP0 * (relativePermittivity[i] - relativePermittivity[0] * sin(theta) * sin(theta)));
		ehNew[i] = ehPrev[i] - prefactorEh * (hhNew[i + 1] - hhNew[i]);

		// Update equation for Ee
		double prefactorEe = timestep / (EP0 * relativePermittivity[i] * gridstep);
		eeNew[i] = eePrev[i] - prefactorEe * (heNew[i+1] - heNew[i]);

		// TFSF correction term
		if (i == tfsfStartIdx)
		{
			double signalTime = timestep * (currTimestep + 0.5) +  pw.Dir_Y  * gridstep / (2 * C0);
			ehNew[i] += prefactorEh * initialUpdateFields[1] * signal_function(signalTime);
			eeNew[i] += prefactorEe * initialUpdateFields[4] * signal_function(signalTime);
		}	
	}
	return {ehNew, eeNew};
}

std::array<std::vector<double>, 2> PlaneWaveOneDimensionalSimulator::update_fields_ycomponents_in_domain(std::vector<double>& ehPrev, std::vector<double>& heNew, std::vector<double>& hyPrev, std::vector<double>& eyPrev)
{
	// vectors for new y-component fields
	std::vector<double> hyNew(gridLength, 0.0);
	std::vector<double> eyNew(gridLength, 0.0);

	for( int i = cpmlNearEndIdx; i < cpmlFarStartIdx + 1; i++)
	{
		hyNew[i] = -hyPrev[i] + 2 * sqrt(EP0 * relativePermittivity[0] / MU0) * sin(theta) * ehPrev[i];
		eyNew[i] = -eyPrev[i] - 2 * sqrt(relativePermittivity[0] * MU0 / EP0) * sin(theta) / relativePermittivity[i] * heNew[i]; 
	}

	return {hyNew, eyNew};
}

/** 
 * @brief Updates the magnetic fields in the CPML region, start and end of domain
 * @param hhNew - vector containing the Hh fields at the next timestep
 * @param heNew - vector containing the He fields at the next timestep
 * @param ehPrev - vector containing the Eh fields at the previous timestep
 * @param hhPrev - vector containing the Hh fields at the previous timestep
 * @param eePrev - vector containing the Ee fields at the previous timestep
 * @param hePrev - vector containing the He fields at the previous timestep
 * @return bool if update was successful
*/
bool PlaneWaveOneDimensionalSimulator::update_magnetic_fields_in_cpml_region(std::vector<double>& hhNew, std::vector<double>& heNew, std::vector<double>& ehPrev, std::vector<double>& hhPrev, std::vector<double>& eePrev, std::vector<double>& hePrev)
{
	// update the new fields in the near CPML region, i.e. before the domain begins
	for (int j = 1; j < gridLength; j++)
	{
		if ( j >= cpmlNearEndIdx && j < cpmlFarStartIdx)
			continue;

		auto [b, c, k] = get_cpml_coefficients_at_index(j);

		double prefactorHe = relativePermittivity[j] * timestep / 
								(MU0 * (relativePermittivity[j] - relativePermittivity[0] * sin(theta) * sin(theta)));

		if (j < cpmlNearEndIdx)
		{
			// update magnetic auxiliary field
			psiHeNear[j] = b * psiHeNear[j] + c * (eePrev[j] - eePrev[j-1]) / gridstep;
			psiHhNear[j] = b * psiHhNear[j] + c * (ehPrev[j] - ehPrev[j-1]) / gridstep;

			// update magnetic field
			heNew[j] = hePrev[j] - prefactorHe * ((eePrev[j] - eePrev[j-1]) / (k * gridstep) + psiHeNear[j]); 
			hhNew[j] = hhPrev[j] - timestep / MU0 * ((ehPrev[j] - ehPrev[j-1]) / (k * gridstep) + psiHhNear[j]);
		}
		if( j >= cpmlFarStartIdx)
		{
			int js = j - cpmlFarStartIdx;

			// update magnetic auxiliary field
			psiHeFar[js] = b * psiHeFar[js] + c * (eePrev[j] - eePrev[j-1]) / gridstep;
			psiHhFar[js] = b * psiHhFar[js] + c * (ehPrev[j] - ehPrev[j-1]) / gridstep;

			// update magnetic field
			heNew[j] = hePrev[j] - prefactorHe * ((eePrev[j] - eePrev[j-1]) / (k * gridstep) + psiHeFar[js]);
			hhNew[j] = hhPrev[j] - timestep / MU0 * ((ehPrev[j] - ehPrev[j-1]) / (k * gridstep) + psiHhFar[js]);
		}

	}
	return true;
}

/**
 * @brief Updates the electric fields in the CPML region, before and after domain
 * @param ehNew - vector containing the Eh fields at the next timestep
 * @param eeNew - vector containing the Ee fields at the next timestep
 * @param hhNew - vector containing the Hh fields at the next timestep
 * @param ehPrev - vector containing the Eh fields at the previous timestep
 * @param heNew - vector containing the He fields at the next timestep
 * @param eePrev - vector containing the Ee fields at the previous timestep
 * @return bool if update was successful
*/
bool PlaneWaveOneDimensionalSimulator::update_electric_fields_in_cpml_region(std::vector<double>& ehNew, std::vector<double>& eeNew ,std::vector<double>& hhNew, std::vector<double>& ehPrev, std::vector<double>& heNew, std::vector<double>& eePrev)
{
	for (int j = 0; j < gridLength - 1; j++)
	{
		if ( j >= cpmlNearEndIdx && j < cpmlFarStartIdx)
			continue;

		auto [b, c, k] = get_cpml_coefficients_at_index(j);

		double prefactorEh = timestep / (EP0 * (relativePermittivity[j] - relativePermittivity[0] * sin(theta) * sin(theta)));

		if (j < cpmlNearEndIdx)
		{
			// update electric auxiliary field
			psiEhNear[j] = b * psiEhNear[j] + c * (hhNew[j+1] - hhNew[j]) / gridstep;
			psiEeNear[j] = b * psiEeNear[j] + c * (heNew[j+1] - heNew[j]) / gridstep;

			// update electric field
			ehNew[j] = ehPrev[j] - prefactorEh * ((hhNew[j+1] - hhNew[j]) / (k * gridstep) + psiEhNear[j]);	
			eeNew[j] = eePrev[j] - timestep / (EP0 * relativePermittivity[j]) * ((heNew[j+1] - heNew[j]) / (k * gridstep) + psiEeNear[j]);
		}

		if (j >= cpmlFarStartIdx)
		{
			int js = j - cpmlFarStartIdx;

			// update electric auxiliary field
			psiEhFar[js] = b * psiEhFar[js] + c * (hhNew[j+1] - hhNew[j]) / gridstep;
			psiEeFar[js] = b * psiEeFar[js] + c * (heNew[j+1] - heNew[j]) / gridstep;

			// update electric field
			ehNew[j] = ehPrev[j] - prefactorEh * ((hhNew[j+1] - hhNew[j]) / (k * gridstep) + psiEhFar[js]);
			eeNew[j] = eePrev[j] -  timestep / (EP0 * relativePermittivity[j]) * ((heNew[j+1] - heNew[j]) / (k * gridstep) + psiEeFar[js]);
		}
	}
	return true;
}

std::array<double, 3> PlaneWaveOneDimensionalSimulator::get_distance_center_to_corner_first_incidence()
{
	std::array<double, 3> dist;

	dist[0] = (pw.Dir_X >= 0) ? _frame.tfsf_n[0] - _frame.focus[0] : _frame.tfsf_f[0] - _frame.focus[0];
	dist[1] = (pw.Dir_Y >= 0) ? _frame.tfsf_n[1] - _frame.focus[1] : _frame.tfsf_f[1] - _frame.focus[1];
	dist[2] = (pw.Dir_Z >= 0) ? _frame.tfsf_n[2] - _frame.focus[2] : _frame.tfsf_f[2] - _frame.focus[2];

	return dist;
}