/**
 * @file materials/hydro/state.cc
 * @brief Definitions of the hydro module state, hydro/state.h
 */
#include "materials/hydro/state.h"

namespace fdtd::materials::hydro
{
using namespace fdtd::materials;
fdtd::array3<double> nd_z_prev;
bool with_pressure;
bool with_mag_field;
bool with_convection;
bool use_backwards_differences;
int nonlocal_barrier_int;
} // namespace fdtd::materials::hydro
