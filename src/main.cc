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
#include "config/errors.h"
#include "grid/state.h"
#include "app/run_state.h"
#include "app/run_config.h"
#include "io/binary.h"
#include "io/text.h"
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
				 "usage: %s COMMAND [-c FILE] [-o DIR] [-i DIR] [-s SECTION.KEY=VALUE]... [-r] [FILE]\n"
				 "\n"
				 "commands:\n"
				 "  simulate   run the simulation\n"
				 "  mesh       build the mesh, write the material maps and the geometry report\n"
				 "  memory     count the array memory of every rank without allocating it\n"
				 "  signal     report the source signal and its spectrum, see signal -h\n"
				 "  check      validate the configuration and print it resolved\n"
				 "  template   print an INI file with every key and its default\n"
				 "\n"
				 "options:\n"
				 "  -c, --config FILE            configuration file\n"
				 "  -o, --output DIR             output directory\n"
				 "  -i, --input DIR              input directory\n"
				 "  -s, --set SECTION.KEY=VALUE  override a configuration key\n"
				 "  -r, --restart                continue an interrupted run from its state files\n"
				 "  -h, --help                   print this help\n"
				 "  -V, --version                print the version\n",
				 program);
}

/// `check`: validate on one process, test the signal against the time steps, print the resolved INI.
int check_configuration(const fdtd::CliOptions &opt)
{
	std::ifstream in(opt.config_file);
	if (!in)
	{
		std::fprintf(stderr, "error: cannot read configuration file '%s'\n", opt.config_file.c_str());
		return 1;
	}
	std::stringstream text;
	text << in.rdbuf();

	auto config = fdtd::config::load_sim_config(fdtd::config::build_config(fdtd::config::with_catalog(text.str(), opt.config_file, opt.overrides), opt.config_file, opt.overrides, opt.input_dir, opt.output_dir));
	fdtd::config::validate_sim_config(config);

	std::printf("%s", fdtd::config::to_ini(config).c_str());
	const std::size_t points = fdtd::config::spectrum_omegas(config.spectrum).size();

	// the signal against the time steps, the same test the run makes before its time loop
	fdtd::config::apply_sim_config(config, 1, true);
	fdtd::app::set_time_step();
	fdtd::app::configure_signal();
	std::string why;
	const int steps = fdtd::app::run_state().nr_max_time_steps;
	const bool fits = fdtd::sources::signal_fits(steps, why);
	if (!fits && !config.source.allow_truncated_signal)
	{
		std::fprintf(stderr, "error: %s\n", why.c_str());
		return 1;
	}
	std::printf("\n; check: ranks,steps,dt_s,dft_points,signal\n");
	std::printf("; check: %d,%d,%.6e,%zu,%s\n", config.parallel.ranks_x * config.parallel.ranks_y * config.parallel.ranks_z, steps, fdtd::app::run_config().time.dt,
				points, fits ? "covered" : "truncated");
	if (!fits)
		std::printf("; warning: %s\n", why.c_str());
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
	const auto &c = driver.config();
	std::vector<unsigned long long> largest(memory::group_count, 0);
	std::vector<unsigned long long> total(memory::group_count, 0);
	unsigned long long bytes_max = 0;
	unsigned long long bytes_min = ~0ull;
	unsigned long long bytes_all = 0;
	for (const auto &rank : per_rank)
	{
		unsigned long long bytes = 0;
		for (int g = 0; g < memory::group_count; ++g)
		{
			bytes += rank[g];
			largest[g] = std::max<unsigned long long>(largest[g], rank[g]);
			total[g] += rank[g];
		}
		bytes_max = std::max(bytes_max, bytes);
		bytes_min = std::min(bytes_min, bytes);
		bytes_all += bytes;
	}
	constexpr double gib = 1073741824.0;
	const int ranks = int(per_rank.size());
	std::printf("layout: ranks_x,ranks_y,ranks_z,ranks\n");
	std::printf("layout: %d,%d,%d,%d\n", c.parallel.ranks_x, c.parallel.ranks_y, c.parallel.ranks_z, ranks);
	std::printf("%s", fdtd::io::memory_table(largest, total, bytes_max, bytes_all).c_str());
	std::printf("job: smallest_rank_gib,overhead_per_rank_gib,largest_rank_gib,total_gib\n");
	std::printf("job: %.3f,%.3f,%.3f,%.3f\n", double(bytes_min) / gib, overhead, double(bytes_max) / gib + overhead, double(bytes_all) / gib + ranks * overhead);
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
		std::fprintf(stderr, "error: unknown command '%s'\n\n", command.c_str());
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
