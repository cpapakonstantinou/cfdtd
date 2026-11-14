#ifndef SOURCES_H
#define SOURCES_H

#include <vector>
#include <string>
#include <array>

//              SOURCE STRUCTS
/// one [dipole.*] section: position in cells and E amplitudes
struct DipoleInput { int x, y, z; double e_x, e_y, e_z; };

struct PlaneWave
{
	double Dir_X;
	double Dir_Y;
	double Dir_Z;
	double AmpEx;
	double AmpEy;
	double AmpEz;
	double AmpHx;
	double AmpHy;
	double AmpHz;
	bool UseDerivative;
	double dispersionFactor = 1.0;
};

struct GaussianBeamData
{
	double maxTheta;
	double maxPhi;
	int nrTheta;
	int nrPhi;
	double deltaTheta;
	double deltaPhi;
	double focalLength;
	double beamWaist;
	double amplConst;
};

/// Where the superposed plane waves are anchored: the focus, the structure centre, the TF/SF box
/// faces, and the delay every wave's signal is shifted by, positions in metres, PML included
struct WaveFrame
{
	double focus[3] = {0, 0, 0};
	double tfsf_n[3] = {0, 0, 0}, tfsf_f[3] = {0, 0, 0}; ///< near and far TF/SF faces
	double offset_delay = 0;                             ///< [s]
};

/// Plane waves of a Gaussian beam and one oblique plane wave, src/sources/gaussian_beam.cc
void build_gaussian_plane_waves(GaussianBeamData &, std::vector<PlaneWave> &);
void build_oblique_plane_wave(double &, double &, double &, bool &,
                           std::vector<PlaneWave> &returnVector);

//              SOURCE INJECTION FUNCTIONS

// Time profile and plane-wave list: sources/signal.h

// Space profile -- plane-wave superposition and Gaussian beam: sources/superposition.h

// dipole sources: sources/source.h

#endif

