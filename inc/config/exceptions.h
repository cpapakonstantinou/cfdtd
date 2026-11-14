#ifndef __EXCEPTIONS_H__
#define __EXCEPTIONS_H__

/** \file   config/exceptions.h
 *  \brief  exception types thrown by the configuration utilities, Config
 *          and by the typed simulation-configuration layer, SimConfig.
 **/
#include <stdexcept>
#include <string>
#include <cstring>

namespace fdtd::config
{
	/// A file could not be opened / read.
	class file_error : public std::runtime_error
	{
		public:
			file_error(const std::string& filename, int err)
			: std::runtime_error("cannot open '" + filename + "': " + std::strerror(err)),
			  _filename(filename), _errno(err) {}

			const std::string& filename() const { return _filename; }
			int error_number() const { return _errno; }

		private:
			std::string _filename;
			int _errno;
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

	/// The configuration is syntactically valid but physically / logically inconsistent.
	class validation_error : public std::runtime_error
	{
		public:
			explicit validation_error(const std::string& msg) : std::runtime_error(msg) {}
	};
} // namespace fdtd::config

#endif // __EXCEPTIONS_H__
