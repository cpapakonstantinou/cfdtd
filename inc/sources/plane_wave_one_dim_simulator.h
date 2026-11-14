#include "sources/signal.h"
#include <vector>
#include <string>
#include <array>
#include "sources/types.h"

class PlaneWaveOneDimensionalSimulator 
{
  public:
	// initialize the simulator for a particular plane wave
	bool init(PlaneWave &, const WaveFrame &frame);

	// run simulation
	bool run_simulation();

	// return fields in cartesian coordinates at given location
	std::array<double, 6> return_fields(double, double, double, double);

	// return medium index at given y position
	int get_medium_index(double);

	// grid properties that need to be public for plotting functionality of OneDimensionalSimulator
	int gridShiftIdx;
	int gridLength;

  private:
	WaveFrame _frame;                    ///< the superposition's anchor, focus, TF/SF faces, delay
	const fdtd::sources::Signal *_signal = nullptr; ///< the run's signal
	// plane wave struct that was used to initialize
	PlaneWave pw;
	// components of the vector defining the principal plane together with y
	double mx, mz;
	// incident angle with y axis
	double theta;
	// tangential phase velocity to media stack
	double tangentialPhaseVelocity;
	// time delay for the individual plane wave
	double timeDelay;

	// grid
	double gridstep;
	double maxPropagationVelocity;
	double courantNumber;
	int numCpmlCells;
	int tfsfStartIdx;
	int cpmlNearEndIdx;
	int cpmlFarStartIdx;
	int numInitialTimesteps;
	int currTimestep;
	double timestep;

	// flag to use derivative of signal function 
	bool useDerivativeSignal;

	// array giving the material index at each grid point
	std::vector<double> relativePermittivity;

	// array giving the medium index at each grid point
	std::vector<int> mediumIndices;

	// initial values in the cartesian and 1D simulation coordinate system 
	std::array<double, 6> initialCartesianFields;
	std::array<double, 6> initialUpdateFields;

	// vectors containing fields at all time steps and locations
	// inner vector for location / outer vector for time steps
	std::vector<std::vector<double>> Eh;
	std::vector<std::vector<double>> Hh;
	std::vector<std::vector<double>> Hy;
	std::vector<std::vector<double>> Ee;
	std::vector<std::vector<double>> He;
	std::vector<std::vector<double>> Ey;
	std::vector<double> times;

	// vectors for CPML auxiliary fields
	std::vector<double> psiEhFar;
	std::vector<double> psiHhFar;
	std::vector<double> psiEeFar;
	std::vector<double> psiHeFar;
	std::vector<double> psiEhNear;
	std::vector<double> psiHhNear;
	std::vector<double> psiEeNear;
	std::vector<double> psiHeNear;

	// setup functions
	bool set_geometry_parameters(PlaneWave&);
	bool set_grid_parameters();
	bool set_velocity_and_timestep();
	bool set_medium_indices();
	bool initialize_field_vectors();

	// transformation functions
	std::array<double, 6>
	cartesian_fields_to_update_fields(const std::array<double, 6>&);
	std::array<double, 6>
	update_fields_to_cartesian_fields(const std::array<double, 6>&);

	// run a single time step
	bool propagate_timestep();
	std::array<std::vector<double>, 2> update_magnetic_fields_in_domain(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
	std::array<std::vector<double>, 2> update_electric_fields_in_domain(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
	std::array<std::vector<double>, 2> update_fields_ycomponents_in_domain(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);

	// calculate coefficients for CPML regions at given location
	std::array<double, 3> get_cpml_coefficients_at_index(int);

	// Update functions for the fields in the CPML region
	bool update_magnetic_fields_in_cpml_region(std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&);
	bool update_electric_fields_in_cpml_region(std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&);

	// signal function 
	double signal_function(double);

	std::array<double, 3> get_distance_center_to_corner_first_incidence();

};
