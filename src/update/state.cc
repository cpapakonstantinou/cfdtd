/**
 * @file update/state.cc
 * @brief Definitions of the update module state, update/state.h
 */
#include "update/state.h"
#include "update/update.h"

namespace fdtd::update
{

} // namespace fdtd::update

namespace fdtd::update
{
namespace
{
Spec the_spec; ///< set once by the application
} // namespace
const Spec &spec() { return the_spec; }
void set_spec(const Spec &s) { the_spec = s; }
} // namespace fdtd::update
