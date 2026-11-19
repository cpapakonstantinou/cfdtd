/** \file   config/errors.h
 *  \brief  The errors of the configuration layer: those of cfg::Config, and validation.
 */
#ifndef FDTD_CONFIG_ERRORS_H
#define FDTD_CONFIG_ERRORS_H

#include <stdexcept>

#include <exceptions.h>

namespace fdtd::config
{
using cfg::file_error;
using cfg::syntax_error;

/// a configuration that parses but is physically or logically inconsistent
class validation_error : public std::runtime_error
{
public:
	using std::runtime_error::runtime_error;
};
} // namespace fdtd::config

#endif // FDTD_CONFIG_ERRORS_H
