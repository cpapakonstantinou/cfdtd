/**
 * @file app/options.h
 * @brief What the command line asked for, filled by main() with getopt
 */
#ifndef FDTD_APP_OPTIONS_H
#define FDTD_APP_OPTIONS_H

#include <string>
#include <vector>

namespace fdtd
{
enum class Command { Simulate, Mesh, Memory, Check, Template, Help, Version };

/// the command and its options
struct CliOptions
{
	Command command = Command::Help;
	std::string config_file;
	std::string output_dir;             ///< empty = keep [paths] output_dir
	std::string input_dir;              ///< empty = keep [paths] input_dir
	std::vector<std::string> overrides; ///< "section.key=value", in command-line order
	bool restart = false;               ///< continue from the state files in <output>/state
};

inline const char *command_name(Command c)
{
	switch (c)
	{
		case Command::Simulate: return "simulate";
		case Command::Mesh:     return "mesh";
		case Command::Memory:   return "memory";
		case Command::Check:    return "check";
		case Command::Template: return "template";
		case Command::Version:  return "version";
		case Command::Help:     break;
	}
	return "help";
}
} // namespace fdtd

#endif
