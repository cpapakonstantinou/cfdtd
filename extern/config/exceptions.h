#ifndef __EXCEPTIONS_H__
#define __EXCEPTIONS_H__

/** \file   exceptions.h
 *  \brief  exception types thrown by the configuration utilities, cfg::Config
 *
 *  Reconstructed: config.h includes this header but it was not among the files received.
 *  It declares what config.cc throws; replace it with the upstream copy.
 **/
#include <stdexcept>
#include <string>
#include <cstring>

namespace cfg
{
	/// A file could not be opened / read.
	class file_error : public std::runtime_error
	{
		public:
			file_error(const std::string& filename, int err)
			: std::runtime_error("cannot open '" + filename + "': " + std::strerror(err)),
			  filename_(filename), errno_(err) {}

			const std::string& filename() const { return filename_; }
			int error_number() const { return errno_; }

		private:
			std::string filename_;
			int errno_;
	};

	/// A configuration file or value is malformed.
	class syntax_error : public std::runtime_error
	{
		public:
			/// Error at a specific line of a file.
			syntax_error(const std::string& filename, unsigned line, const std::string& msg = "syntax error")
			: std::runtime_error(filename + ":" + std::to_string(line) + ": " + msg) {}

			/// Error on a specific [section] key.
			syntax_error(const std::string& sect, const std::string& key, const std::string& msg)
			: std::runtime_error("[" + sect + "] " + key + ": " + msg) {}
	};

} // namespace cfg

#endif // __EXCEPTIONS_H__
