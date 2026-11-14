#include "core/constants.h"
#include "materials/state.h"
#include "sources/types.h"
#include "io/binary.h"

#include <string>
#include <vector>
#include "sources/spec.h"
using fdtd::constants::INVVACIMP;
using fdtd::constants::PI2;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
} // namespace

/**
 * @brief Iterate over all angles and construct plane waves
 * 
 * @param returnVector [out] Address of vector to store PlaneWaves in.
 */
void build_gaussian_plane_waves(GaussianBeamData &beamData, std::vector<PlaneWave> &returnVector)
{
	int maxPhiLocal;
	bool useDerivative = true;

	for (int theta_index = 0; theta_index <= beamData.nrTheta - 1; ++theta_index)
	{
		// determine current theta and maximum index for phi
		double theta = theta_index * beamData.deltaTheta;
		maxPhiLocal = (theta_index == 0) ? 0 : beamData.nrPhi - 1;

		for (int phi_index = 0; phi_index <= maxPhiLocal; ++phi_index)
		{
			// determine current phi
			double AMPL;
			double phi = phi_index * beamData.deltaPhi;

			if (sp.cw_pulse)
				AMPL = 1 / double(PI2 * sp.highest_frequency);
			else
			{
				double exponent = beamData.focalLength * sin(theta) / beamData.beamWaist;
				AMPL = beamData.amplConst * sqrt(cos(theta)) * sin(theta) * exp(-exponent * exponent);
			}

			build_oblique_plane_wave(theta, phi, AMPL, useDerivative, returnVector);
			
		}
	}
}

/**
 * @brief build single PlaneWave instance for a given theta, phi, and amplitude
 * 
 * @param theta [in] Angle with respect to y axis
 * @param phi [in] Angle with respect to x axis in the x-z plane
 * @param AMP [in] Amplitude of the plane wave instance
 * @param returnVector [out] Address of the vector of plane waves
 */
void build_oblique_plane_wave(double &theta, double &phi, double &AMPL, bool &useDerivative, std::vector<PlaneWave> &returnVector)
{
	const auto &mp = fdtd::materials::model_parameters();
	// push back vector to add new plane wave
	returnVector.push_back(PlaneWave());

	// reduce function calls by pre-assigning these
	double sin_theta = sin(theta);
	double cos_theta = cos(theta);
	double sin_phi = sin(phi);
	double cos_phi = cos(phi);

	// address so that the following lines are shorter
	PlaneWave &tmpPlaneWave = returnVector[returnVector.size() - 1];

	// Determine plane wave direction vector
	tmpPlaneWave.Dir_X = -sin_theta * cos_phi;
	tmpPlaneWave.Dir_Y = cos_theta;
	tmpPlaneWave.Dir_Z = -sin_theta * sin_phi;
	
	// Calculate plane wave amplitudes
	tmpPlaneWave.AmpEx = AMPL * sin_phi * cos_phi * (cos_theta - 1.);
	tmpPlaneWave.AmpEy = AMPL * sin_theta * sin_phi;
	tmpPlaneWave.AmpEz = AMPL * (cos_phi * cos_phi + cos_theta * sin_phi * sin_phi);

	double scaleH = INVVACIMP * mp.n_refr_index[sp.medium[0]];
	tmpPlaneWave.AmpHx = (tmpPlaneWave.Dir_Y * tmpPlaneWave.AmpEz
						  - tmpPlaneWave.Dir_Z * tmpPlaneWave.AmpEy) * scaleH;
	tmpPlaneWave.AmpHy = (tmpPlaneWave.Dir_Z * tmpPlaneWave.AmpEx
						  - tmpPlaneWave.Dir_X * tmpPlaneWave.AmpEz) * scaleH;
	tmpPlaneWave.AmpHz = (tmpPlaneWave.Dir_X * tmpPlaneWave.AmpEy
						  - tmpPlaneWave.Dir_Y * tmpPlaneWave.AmpEx) * scaleH;

	tmpPlaneWave.UseDerivative = useDerivative;
	return;
}
