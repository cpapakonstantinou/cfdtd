#ifndef __APP_MACHINE_H__
#define __APP_MACHINE_H__
/** \file   app/machine.h
 *  \brief  Application-level state machine, boost-ext/sml, of the 3D-FDTD driver.
 *
 *  The machine is a template over the driver type D so that the control flow
 *  can be unit-tested with a fake driver, tests/test_machine.cpp. D provides
 *  the phase actions, each one collective over MPI, recording success via
 *  ok(), and the guards used to branch.
 *
 *            +-----------------------------------------------------------+
 *            |                         [failed] --> report --> X         |
 *            v                                                           |
 *   Boot --> Configured --> Decomposed --> Initialized --> Ready ---+
 *                                                           |  |  |
 *                         [mesh_only] dump_geometry  <------+  |  |
 *                         [dummy_run] write_outputs  <---------+  |
 *                                       time_loop    <------------+
 *                                           |
 *                                        Stepped --finish_run--> Done
 *                                                                   |   |
 *                                                    write_outputs <----+
 *                                                         |
 *                                                     Finishing --finalize--> X
 *
 *  Every transition is triggered by the same `next` event: the driver loop just
 *  feeds `next` until the machine terminates. What happens is decided by the
 *  current state and the guards, never by flags threaded through main().
 **/
#include <boost/sml.hpp>

namespace fdtd
{
namespace sm
{
	/// The only event: "the current phase is done, advance".
	struct next {};

	// States
	struct Boot {};
	struct Configured {};
	struct Decomposed {};
	struct Initialized {};
	struct Ready {};
	struct Restored {}; ///< the state files of an interrupted run are loaded
	struct Stepped {};
	struct Done {};
	struct Interrupted {}; ///< a signal stopped the time loop: write the outputs it has
	struct Finishing {};

	template <class D>
	struct AppMachine
	{
		auto operator()() const
		{
			using namespace boost::sml;

			// guards
			const auto failed      = [](const D& d) { return !d.ok(); };
			const auto mesh_only   = [](const D& d) { return d.mesh_only(); };
			const auto dummy_run   = [](const D& d) { return d.dummy_run(); };
			const auto interrupted = [](const D& d) { return d.interrupted() != 0; };
			const auto restart     = [](const D& d) { return d.restart(); };

			// actions
			const auto load_config   = [](D& d) { d.load_configuration(); };
			const auto decompose     = [](D& d) { d.decompose(); };
			const auto init_ranks    = [](D& d) { d.initialize_ranks(); };
			const auto prepare_run = [](D& d) { d.prepare_run(); };
			const auto dump_geometry = [](D& d) { d.dump_geometry(); };
			const auto load_state    = [](D& d) { d.load_state(); };
			const auto time_loop     = [](D& d) { d.run_time_loop(); };
			const auto finish_run  = [](D& d) { d.finish_run(); };
			const auto write_outputs = [](D& d) { d.write_outputs(); };
			const auto finalize      = [](D& d) { d.finalize(); };
			const auto report        = [](D& d) { d.report_error(); };

			// clang-format off
			return make_transition_table(
				*state<Boot>        + event<next>                 / load_config   = state<Configured>,

				 state<Configured>  + event<next> [failed]        / report        = X,
				 state<Configured>  + event<next>                 / decompose     = state<Decomposed>,

				 state<Decomposed>  + event<next> [failed]        / report        = X,
				 state<Decomposed>  + event<next>                 / init_ranks    = state<Initialized>,

				 state<Initialized> + event<next> [failed]        / report        = X,
				 state<Initialized> + event<next>                 / prepare_run = state<Ready>,

				 state<Ready>  + event<next> [failed]        / report        = X,
				 state<Ready>  + event<next> [mesh_only]     / dump_geometry = state<Finishing>,
				 state<Ready>  + event<next> [dummy_run]     / write_outputs = state<Finishing>,
				 state<Ready>  + event<next> [restart]       / load_state    = state<Restored>,
				 state<Ready>  + event<next>                 / time_loop     = state<Stepped>,

				 state<Restored>    + event<next> [failed]        / report        = X,
				 state<Restored>    + event<next>                 / time_loop     = state<Stepped>,

				 state<Stepped>     + event<next> [failed]        / report        = X,
				 state<Stepped>     + event<next> [interrupted]   / finish_run  = state<Interrupted>,
				 state<Stepped>     + event<next>                 / finish_run  = state<Done>,

				 state<Interrupted> + event<next> [failed]        / report        = X,
				 state<Interrupted> + event<next>                 / write_outputs = state<Finishing>,

				 state<Done>   + event<next> [failed]        / report        = X,
				 state<Done>   + event<next>                 / write_outputs = state<Finishing>,

				 state<Finishing>   + event<next> [failed]        / report        = X,
				 state<Finishing>   + event<next>                 / finalize      = X
			);
			// clang-format on
		}
	};

	/** \brief Drive the machine to completion. \returns true when the run ended without error. */
	template <class D>
	bool run_machine(D& driver)
	{
		boost::sml::sm<AppMachine<D>> machine{driver};
		while (!machine.is(boost::sml::X))
			machine.process_event(next{});
		return driver.ok();
	}
} // namespace sm
} // namespace fdtd

#endif // __APP_MACHINE_H__
