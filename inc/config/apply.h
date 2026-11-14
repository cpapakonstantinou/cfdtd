#ifndef __APPLY_CONFIG_H__
#define __APPLY_CONFIG_H__
/** \file   config/apply.h
 *  \brief  Mapping of the typed configuration onto the legacy solver state.
 **/
#include "config/sim_config.h"

namespace fdtd::config
{
/** \brief Fill app::RunConfig and the remaining solver globals from the configuration.
 *
 *  Must be called on EVERY rank with the same configuration, replaces the
 *  root-only ReadJSONDataFile / ReadIniDataFile.
 *  \throws validation_error, e.g. unreadable chi tensor file */
void apply_sim_config(const SimConfig& c, int world_size, bool is_root);
} // namespace fdtd::config

#endif // __APPLY_CONFIG_H__
