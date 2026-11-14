/**
 * @file io/text.h
 * @brief Text output and input helpers: the run log and the legacy table reader
 */
#ifndef FDTD_IO_TEXT_H
#define FDTD_IO_TEXT_H

#include <ctime>
#include <fstream>
#include <string>

namespace fdtd::io
{
/// appends text to <dirPath>/<logName><logNumber>.txt, with the time since startTime; also printed when verbose
bool write_log(std::string text, std::string dirPath, clock_t startTime, std::string logName = "log", bool verboseFlag = false,
			   int logNumber = 0);
/**
 * @brief One log record per line on stdout: `key: text`, or `key[rank]: text` from a rank other than 0
 *
 * The key names the module or the quantity, so `grep '^step:' log` extracts a table. A message with
 * several lines gives one record per line; blank lines and lines of dashes are dropped.
 */
void log(const char *key, const char *format, ...) __attribute__((format(printf, 2, 3)));
inline bool log_muted = false; ///< log() prints nothing, `fdtd memory`: the set-up runs once per rank
/// the first space-separated word of the next line
std::string read_text_entry(std::fstream &file);
} // namespace fdtd::io

#endif
