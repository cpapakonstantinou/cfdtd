/**
 * @file materials/nonlinear/state.cc
 * @brief Definitions of the nonlinear module state, nonlinear/state.h
 */
#include "materials/nonlinear/state.h"

namespace fdtd::materials::nonlinear
{
using namespace fdtd::materials;

double chi_2[max_slots][18];
double chi_3[max_slots][21];
bool iso_chi_2[max_slots];
bool iso_chi_3[max_slots];
bool nl[max_slots];
} // namespace fdtd::materials::nonlinear
