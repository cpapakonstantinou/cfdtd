#include <vector>
#include <string>
#include <array>
#include "sources/types.h"
#include "sources/plane_wave_one_dim_simulator.h"

class OneDimensionalSimulator
{
	public:
		// initialize the simulator for a vector fo plane waves
		bool init(std::vector<PlaneWave>, const WaveFrame &frame);
		const WaveFrame &frame() const { return _frame; } ///< where the waves are anchored

		// write the results of the 1D simulations
		bool write_onedim_results();

		// return the fields in cartesian coordinates at given time and location
		std::array<double, 6> return_fields(double, double, double, double);

		// returns the medium index at a given y location
		int get_medium_index(double);

	private:
		// list of simulators for single plane waves
		std::vector<PlaneWaveOneDimensionalSimulator> singleSimulatorList;
		WaveFrame _frame;

		void check_media_ordering();
};

