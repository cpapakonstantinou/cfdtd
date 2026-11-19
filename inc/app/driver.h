#ifndef __DRIVER_H__
#define __DRIVER_H__
/** \file   app/driver.h
 *  \brief  The 3D-FDTD application driver: one method per phase of a run.
 *
 *  The phases are sequenced by the state machine in app_machine.h; this class
 *  only knows how to perform each phase and to report whether it succeeded.
 *  Every phase is collective: all ranks call it, and the success flag is
 *  agreed on with an MPI reduction, so all ranks always take the same
 *  transition.
 **/
#include <memory>
#include <string>

#include "app/options.h"
#include "core/memory.h"
#include <array>
#include <vector>
#include "config/sim_config.h"

namespace fdtd
{

class Driver
{
	public:
		Driver(const CliOptions& opt, int rank, int size);
		~Driver();
		Driver(const Driver&) = delete;
		Driver& operator=(const Driver&) = delete;

		// guards, see app_machine.h
		bool ok() const { return _ok; }
		bool mesh_only() const { return _opt.command == Command::Mesh || _opt.command == Command::Memory; }
		bool dummy_run() const { return _dummy_run; }
		bool restart() const { return _opt.restart; }
		int interrupted() const { return _signal; } ///< signal that stopped the time loop, all ranks agree, or 0

		// phases
		void load_configuration(); ///< read + broadcast INI, apply CLI overrides, validate, fill the lattice parameters and globals, create output dirs
		void decompose();         ///< root: domain setup + MPI decomposition; then broadcast to all ranks
		void initialize_ranks();   ///< per-rank subdomain initialization
		void prepare_run();      ///< pulse, allocation, geometry, sources and monitors of the run
		void dump_geometry();      ///< mesh path: material index dumps + geometry report
		void load_state();         ///< restart: load the state files of an interrupted run
		void run_time_loop();       ///< FDTD iterations of the run
		void finish_run();       ///< DFT normalization and the analysis, ART, Q, exact multipoles, ...
		void write_outputs();      ///< bulk DFT output, custom monitors, release memory
		void finalize();          ///< final log line, optimizer shutdown
		void report_error();       ///< print the agreed failure

		/// `fdtd memory`: bytes of the arrays each rank of the configured layout would allocate before
		/// the geometry is built, per group, core::memory::Group; one process plays every rank in turn.
		bool estimate_memory(std::vector<std::array<unsigned long long, core::memory::group_count>>& per_rank);

		const config::SimConfig& config() const { return _cfg; }

	private:
		/// Collective success check: false on any rank -> false on all ranks.
		bool agree(bool local_ok, const std::string& what);

		// helpers of PrepareRun()
		void configure_run();
		bool allocate_sweep_memory();
		void build_geometry();
		void setup_sources_and_monitors();
		void build_run_lists(); ///< the source and monitor lists of the run
		void write_geometry_report();

		// helpers of run_time_loop()
		void update_magnetic();
		void update_electric();
		void accumulate_monitors();
		void print_trace();
		bool wallclock_exceeded();
		void write_time_domain_frames();
		void periodic_analysis();

		struct Context; ///< legacy-typed per-run state, defined in driver.cpp
		std::unique_ptr<Context> _ctx;

		CliOptions _opt;
		config::SimConfig _cfg;
		int _rank;
		int _size;
		bool _ok = true;
		bool _dummy_run = false;
		int _signal = 0; ///< see interrupted()
		bool _state_written = false; ///< this run ended with a state dump
		int _completed_steps = 0; ///< steps done before this time loop, restart
		void save_state();         ///< dump the state after an interrupted time loop
		std::string _error;
};

} // namespace fdtd

#endif // __DRIVER_H__
