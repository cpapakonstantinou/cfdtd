/** \file io/report.cc
 *  \brief Text reports: the run's metadata as comment lines above the column header. */
#include "io/text.h"
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
