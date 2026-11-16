// Unit tests for the application state machine with a fake driver, no MPI, no solver.
#include "app/machine.h"

#include <cstdio>
#include <string>
#include <vector>

static int failures = 0;
#define CHECK(cond)                                                              \
	do {                                                                         \
		if (!(cond)) { std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); ++failures; } \
	} while (0)

struct FakeDriver
{
	std::vector<std::string> log;
	std::string fail_in;      // name of the phase that fails
	bool mesh = false;
	bool dummy = false;
	int sweep = 0;
	bool _ok = true;

	void step(const std::string& name) { log.push_back(name); if (name == fail_in) _ok = false; }

	bool ok() const { return _ok; }
	bool mesh_only() const { return mesh; }
	bool dummy_run() const { return dummy; }
	int interrupt_in_sweep = -1; ///< the time loop of this sweep is stopped by a signal
	int signal = 0;
	int interrupted() const { return signal; }
	bool from_state = false;
	bool restart() const { return from_state; }
	void load_state() { step("restore"); }

	void load_configuration() { step("load"); }
	void decompose() { step("decompose"); }
	void initialize_ranks() { step("init"); }
	void prepare_run() { step("prepare" + std::to_string(sweep)); }
	void dump_geometry() { step("dump"); }
	void run_time_loop() { step("loop" + std::to_string(sweep)); if (sweep == interrupt_in_sweep) signal = 15; }
	void finish_run() { step("finish" + std::to_string(sweep)); }
	void write_outputs() { step("write"); }
	void finalize() { step("finalize"); }
	void report_error() { log.push_back("error"); }
};

static std::string join(const std::vector<std::string>& v)
{
	std::string s;
	for (const auto& x : v) s += (s.empty() ? "" : " ") + x;
	return s;
}

int main()
{
	{ // simulate
		FakeDriver d;
		CHECK(fdtd::sm::run_machine(d));
		CHECK(join(d.log) == "load decompose init prepare0 loop0 finish0 write finalize");
	}
	{ // a signal in the time loop: finish the run, write, finalize
		FakeDriver d;
		d.interrupt_in_sweep = 0;
		CHECK(fdtd::sm::run_machine(d));
		CHECK(join(d.log) == "load decompose init prepare0 loop0 finish0 write finalize");
	}
	{ // restart: the state is loaded between the preparation and the time loop
		FakeDriver d;
		d.from_state = true;
		CHECK(fdtd::sm::run_machine(d));
		CHECK(join(d.log) == "load decompose init prepare0 restore loop0 finish0 write finalize");
		FakeDriver bad;
		bad.from_state = true;
		bad.fail_in = "restore";
		CHECK(!fdtd::sm::run_machine(bad) && bad.log.back() == "error");
	}
	{ // mesh path never runs the time loop
		FakeDriver d;
		d.mesh = true;
		CHECK(fdtd::sm::run_machine(d));
		CHECK(join(d.log) == "load decompose init prepare0 dump finalize");
	}
	{ // mode-source dummy run skips straight to output
		FakeDriver d;
		d.dummy = true;
		CHECK(fdtd::sm::run_machine(d));
		CHECK(join(d.log) == "load decompose init prepare0 write finalize");
	}
	{ // failures stop the machine right after the failing phase
		for (const char* phase : {"load", "decompose", "init", "prepare0", "loop0", "finish0", "write"})
		{
			FakeDriver d;
			d.fail_in = phase;
			CHECK(!fdtd::sm::run_machine(d));
			CHECK(d.log.back() == "error");
			CHECK(d.log[d.log.size() - 2] == phase);
		}
	}
	std::printf(failures ? "%d FAILURE(S)\n" : "all state-machine tests passed\n", failures);
	return failures ? 1 : 0;
}
