/**
 * @file io/text.h
 * @brief Text output and input helpers: the run log and the legacy table reader
 */
#ifndef FDTD_IO_TEXT_H
#define FDTD_IO_TEXT_H

#include <cstdio>
#include <ctime>
#include <fstream>
#include <string>
#include <vector>

namespace fdtd::io
{
/// appends text to <dir>/run.log, with the time since start; also a `run:` record when verbose
bool write_log(const std::string &text, const std::string &dir, clock_t start, bool verbose = false);
/// the lines, without "# ", written above the column header of every report: what the run was
void set_report_metadata(std::vector<std::string> lines);
/**
 * @brief Opens a text report for writing: the metadata, then the column header, each line after "# "
 *
 * Readers skip the comment lines: `numpy.loadtxt(path)`, `pandas.read_csv(path, comment="#", sep=r"\s+", header=None)`.
 * @return the open file, nullptr when it cannot be opened
 */
FILE *open_report(const std::string &path, const std::string &columns);
/**
 * @brief One log record per line on stdout: `key: text`, or `key[rank]: text` from a rank other than 0
 *
 * The key names the module or the quantity, so `grep '^step:' log` extracts a table. A message with
 * several lines gives one record per line; blank lines and lines of dashes are dropped.
 */
void log(const char *key, const char *format, ...) __attribute__((format(printf, 2, 3)));
/**
 * @brief The array memory as two tables: `memory:` per group and in all, largest rank and total in
 *        GiB, then `bytes:` with the exact byte counts of all arrays
 *
 * `fdtd memory` and the geometry report of `fdtd mesh` print the same tables, so they compare line by line.
 */
std::string memory_table(const std::vector<unsigned long long> &largest, const std::vector<unsigned long long> &total, unsigned long long largest_all,
						 unsigned long long total_all);
inline bool log_muted = false; ///< log() prints nothing, `fdtd memory`: the set-up runs once per rank
/// the first space-separated word of the next line
std::string read_text_entry(std::fstream &file);
} // namespace fdtd::io

#endif
