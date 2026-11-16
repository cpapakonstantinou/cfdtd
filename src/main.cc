/** \file   main.cc
 *  \brief  Entry point of 3D-FDTD: the command line, getopt, and MPI.
 *
 *    fdtd template                      > run.ini     # every key, documented, with defaults
 *    fdtd check    run.ini [-s ...]                   # validate, test the signal, print the resolved config, no MPI
 *    fdtd memory run.ini                              # memory per rank and in total, one process
 *    mpirun -n N fdtd mesh     run.ini -o out/        # build + dump the mesh, geometry report
 *    mpirun -n N fdtd simulate run.ini -o out/        # full run
 *
 *  All sequencing lives in the state machine, app/machine.h; all phases in fdtd::Driver, app/driver.cc.
 **/
#include "app/signal_report.h"
#include "app/options.h"
#include "app/setup.h"
#include "app/signals.h"
#include "app/machine.h"
#include "config/apply.h"
#include "app/driver.h"
#include "config/exceptions.h"
#include "grid/state.h"
#include "app/run_state.h"
#include "app/run_config.h"
#include "io/binary.h"
#include "sources/one_dim_simulator.h"
#include "sources/signal.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <unistd.h>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include "app/signal_config.h"
// no shared state

#ifndef FDTD_VERSION
#define FDTD_VERSION "dev"
#endif

namespace
{
void print_usage(std::FILE *to, const char *program)
{
	std::fprintf(to,
				 "usage: %s <command> [options] [config.ini]\n"
				 "\n"
				 "commands:\n"
				 "  simulate   build the mesh and run the FDTD simulation\n"
				 "  mesh       build the mesh, write the material index maps and the geometry report\n"
				 "  memory     count the memory the run needs, per rank and in total, without allocating it;\n"
				 "             one process, no mpirun needed\n"
				 "  signal     report the signal length, the time steps it needs and its spectrum, no MPI;\n"
				 "             see fdtd signal --help\n"
				 "  check      validate the configuration, test that the time steps cover the signal,\n"
				 "             print the resolved INI, no MPI\n"
				 "  template   print a commented INI template with every key and its default\n"
				 "\n"
				 "options:\n"
				 "  -c, --config FILE          configuration file, or give it after the command\n"
				 "  -o, --output DIR           output directory, overrides [paths] output_dir\n"
				 "  -i, --input DIR            input/data directory, overrides [paths] input_dir\n"
				 "  -s, --set SECTION.KEY=VAL  override a configuration key, repeatable\n"
				 "  -r, --restart              simulate: continue from the state files an interrupted run left\n"
				 "                             in <output>/state, same configuration and rank layout\n"
				 "  -h, --help                 show this help\n"
				 "  -V, --version              show the version\n"
				 "\n"
				 "precedence: config file < --set < -o/-i\n"
				 "example:\n"
				 "  mpirun -n 8 %s simulate run.ini -o /scratch/run42 -s grid.dx=2e-9 -s spectrum.points=51\n",
				 program, program);
}

/// `check`: validate on one process, test the signal against the time steps, print the resolved INI.
int check_configuration(const fdtd::CliOptions &opt)
{
	std::ifstream in(opt.config_file);
	if (!in)
	{
		std::fprintf(stderr, "cannot read configuration file '%s'\n", opt.config_file.c_str());
		return 1;
	}
	std::stringstream text;
	text << in.rdbuf();

	auto config = fdtd::config::load_sim_config(fdtd::config::build_config(text.str(), opt.config_file, opt.overrides, opt.input_dir, opt.output_dir));
	fdtd::config::validate_sim_config(config);

	std::printf("%s", fdtd::config::to_ini(config).c_str());
	const auto omegas = fdtd::config::spectrum_omegas(config.spectrum);
	std::printf("\n; spectrum: %zu DFT points\n", omegas.size());
	for (double w : omegas)
		std::printf(";   %.6f nm\n", 2.0 * 3.141592653589793 * 299792458.0 / w * 1e9);

	// the signal against the time steps, the same test the run makes before its time loop
	fdtd::config::apply_sim_config(config, 1, true);
	fdtd::app::set_time_step();
	fdtd::app::configure_signal();
	std::string why;
	const int steps = fdtd::app::run_state().nr_max_time_steps;
	if (!fdtd::sources::signal_fits(steps, why))
	{
		if (!config.source.allow_truncated_signal)
		{
			std::fprintf(stderr, "error: %s\n", why.c_str());
			return 1;
		}
		std::printf("; warning: %s\n", why.c_str());
	}
	else
		std::printf("; signal: covered by the %d time steps (dt = %.6e s)\n", steps, fdtd::app::run_config().time.dt);
	std::printf("; configuration is valid (MPI ranks required: %d)\n",
				config.parallel.ranks_x * config.parallel.ranks_y * config.parallel.ranks_z);
	return 0;
}
/// Resident memory of this process now [MiB], 0 if /proc is not available.
double resident_mib()
{
	long pages = 0;
	long resident = 0;
	if (std::FILE *f = std::fopen("/proc/self/statm", "r"))
	{
		if (std::fscanf(f, "%ld %ld", &pages, &resident) != 2)
			resident = 0;
		std::fclose(f);
	}
	return double(resident) * double(sysconf(_SC_PAGESIZE)) / 1048576.0;
}

/// `memory`: the arrays every rank of the configured layout allocates, counted without allocating them.
int estimate_memory(const fdtd::CliOptions &opt)
{
	namespace memory = fdtd::core::memory;
	// What a rank holds besides its arrays, MPI, HDF5, the executable, small containers: 40 to 65 MiB
	// measured in real runs with Open MPI; a process started without mpirun, as this one, shows less.
	const double overhead = std::max(60.0, resident_mib()) / 1024.0; // [GiB]
	fdtd::Driver driver(opt, 0, 1);
	std::vector<std::array<unsigned long long, memory::group_count>> per_rank;
	if (!driver.estimate_memory(per_rank))
	{
		driver.report_error();
		return 1;
	}
	constexpr double GIB = 1073741824.0;
	const auto &c = driver.config();
	std::array<double, memory::group_count> largest{};
	std::array<double, memory::group_count> total{};
	double rank_max = 0;
	double rank_min = 1e300;
	double all = 0;
	unsigned long long bytes_max = 0;
	unsigned long long bytes_all = 0;
	for (const auto &rank : per_rank)
	{
		unsigned long long bytes = 0;
		for (unsigned long long b : rank)
			bytes += b;
		bytes_max = std::max(bytes_max, bytes), bytes_all += bytes;
		double sum = 0;
		for (int g = 0; g < memory::group_count; ++g)
		{
			const double gib = double(rank[g]) / GIB;
			largest[g] = std::max(largest[g], gib), total[g] += gib, sum += gib;
		}
		rank_max = std::max(rank_max, sum), rank_min = std::min(rank_min, sum), all += sum;
	}
	const int ranks = int(per_rank.size());
	std::printf("memory of %s: %d x %d x %d = %d ranks\n\n", opt.config_file.c_str(), c.parallel.ranks_x, c.parallel.ranks_y, c.parallel.ranks_z, ranks);
	std::printf("arrays [GiB]             largest rank           total\n");
	for (int g = 0; g < memory::group_count; ++g)
		std::printf("  %-18s %16.3f %15.3f\n", memory::GROUP_NAMES[g], largest[g], total[g]);
	std::printf("  %-18s %16.3f %15.3f\n", "all arrays", rank_max, all);
	std::printf("  smallest rank      %16.3f\n", rank_min);
	std::printf("  in bytes           %16llu %15llu\n\n", bytes_max, bytes_all);
	std::printf("process overhead   : about %.2f GiB per rank (MPI, HDF5 and other libraries; 0.04 to 0.065 GiB measured with Open MPI)\n", overhead);
	std::printf("job                : about %.2f GiB on the largest rank, %.2f GiB in total\n\n", rank_max + overhead, all + ranks * overhead);
	std::printf("The arrays are counted by running the solver's own set-up for each rank, up to the mesh, with\n"
				"arrays that take no memory until written: the figures equal those of the geometry report of\n"
				"`fdtd mesh`. Not included: buffers a run allocates after the mesh (a few percent at most in\n"
				"the cases measured; the pre-computed fields of layered plane-wave sources can be larger).\n"
				"Try another layout with -s parallel.ranks_x=.. -s parallel.ranks_y=.. -s parallel.ranks_z=..\n");
	return 0;
}
} // namespace

int main(int argc, char *argv[])
{
	const char *program = argc > 0 ? argv[0] : "fdtd";
	if (argc > 1 && std::string(argv[1]) == "signal")
		return fdtd::app::signal_report(argc - 1, argv + 1);

	// command line
	fdtd::CliOptions opt;
	static const option long_options[] = {{"config", required_argument, nullptr, 'c'}, {"output", required_argument, nullptr, 'o'},
										  {"input", required_argument, nullptr, 'i'},  {"set", required_argument, nullptr, 's'},
										  {"restart", no_argument, nullptr, 'r'},
										  {"help", no_argument, nullptr, 'h'},         {"version", no_argument, nullptr, 'V'},
										  {nullptr, 0, nullptr, 0}};
	for (int c; (c = getopt_long(argc, argv, "c:o:i:s:rhV", long_options, nullptr)) != -1;)
		switch (c)
		{
			case 'c': opt.config_file = optarg; break;
			case 'o': opt.output_dir = optarg; break;
			case 'i': opt.input_dir = optarg; break;
			case 's': opt.overrides.push_back(optarg); break;
			case 'r': opt.restart = true; break;
			case 'h': print_usage(stdout, program); return 0;
			case 'V': std::printf("3D-FDTD %s\n", FDTD_VERSION); return 0;
			default: print_usage(stderr, program); return 2; // getopt printed the reason
		}
	// what is left: the command, then the configuration file
	if (optind >= argc)
	{
		print_usage(stdout, program);
		return 0;
	}
	const std::string command = argv[optind++];
	if (command == "simulate" || command == "sim" || command == "run")
		opt.command = fdtd::Command::Simulate;
	else if (command == "mesh")
		opt.command = fdtd::Command::Mesh;
	else if (command == "memory")
		opt.command = fdtd::Command::Memory;
	else if (command == "check")
		opt.command = fdtd::Command::Check;
	else if (command == "template")
		opt.command = fdtd::Command::Template;
	else if (command == "help")
		opt.command = fdtd::Command::Help;
	else
	{
		std::fprintf(stderr, "error: unknown command '%s' (expected simulate, mesh, memory, check or template)\n\n", command.c_str());
		print_usage(stderr, program);
		return 2;
	}
	if (optind < argc && opt.config_file.empty())
		opt.config_file = argv[optind++];
	if (optind < argc)
	{
		std::fprintf(stderr, "error: unexpected argument '%s'\n\n", argv[optind]);
		print_usage(stderr, program);
		return 2;
	}

	// commands that never touch MPI
	if (opt.command == fdtd::Command::Help)
	{
		print_usage(stdout, program);
		return 0;
	}
	if (opt.command == fdtd::Command::Template)
	{
		std::printf("%s", fdtd::config::template_ini().c_str());
		return 0;
	}
	if (opt.config_file.empty())
	{
		std::fprintf(stderr, "error: '%s' needs a configuration file\n\n", fdtd::command_name(opt.command));
		print_usage(stderr, program);
		return 2;
	}
	if (opt.command == fdtd::Command::Check)
	{
		try
		{
			return check_configuration(opt);
		}
		catch (const std::exception &e)
		{
			std::fprintf(stderr, "error: %s\n", e.what());
			return 1;
		}
	}

	// mesh / simulate
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &fdtd::grid::world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &fdtd::grid::world_rank);

	if (opt.command == fdtd::Command::Memory)
	{
		const int code = estimate_memory(opt);
		MPI_Finalize();
		return code;
	}

	fdtd::app::install_signal_handlers();
	int code = 1;
	{
		fdtd::Driver driver(opt, fdtd::grid::world_rank, fdtd::grid::world_size);
		const bool success = fdtd::sm::run_machine(driver);
		code = !success ? 1 : driver.interrupted() ? 128 + driver.interrupted() : 0;
	}

	MPI_Finalize();
	return code;
}
