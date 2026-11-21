/** \file io/report.cc
 *  \brief Text reports: the run's metadata as comment lines above the column header. */
#include "io/text.h"
#include "core/memory.h"

#include <cstdio>
#include <utility>

namespace
{
std::vector<std::string> metadata;
} // namespace

void fdtd::io::set_report_metadata(std::vector<std::string> lines)
{
	metadata = std::move(lines);
}

FILE *fdtd::io::open_report(const std::string &path, const std::string &columns)
{
	FILE *f = std::fopen(path.c_str(), "w");
	if (!f)
		return nullptr;
	for (const std::string &line : metadata)
		std::fprintf(f, "# %s\n", line.c_str());
	if (!columns.empty())
		std::fprintf(f, "# %s\n", columns.c_str());
	return f;
}

std::string fdtd::io::memory_table(const std::vector<unsigned long long> &largest, const std::vector<unsigned long long> &total, unsigned long long largest_all,
								   unsigned long long total_all)
{
	constexpr double gib = 1073741824.0;
	char line[160];
	std::string t = "memory: group,largest_rank_gib,total_gib\n";
	for (int g = 0; g < fdtd::core::memory::group_count; ++g)
	{
		std::snprintf(line, sizeof line, "memory: %s,%.3f,%.3f\n", fdtd::core::memory::GROUP_NAMES[g], double(largest[g]) / gib, double(total[g]) / gib);
		t += line;
	}
	std::snprintf(line, sizeof line, "memory: all,%.3f,%.3f\nbytes: largest_rank,total\nbytes: %llu,%llu\n", double(largest_all) / gib, double(total_all) / gib,
				  largest_all, total_all);
	return t + line;
}
