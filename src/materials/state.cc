/**
 * @file materials/state.cc
 * @brief Definitions of the materials module state, materials/state.h
 */
#include "materials/state.h"
#include "materials/spec.h"
#include "materials/dispersive_state.h"

namespace fdtd::materials
{
bool laplacian;
bool gradient_divergence;
fdtd::array3<double> magnetization_x;
fdtd::array3<double> magnetization_y;
fdtd::array3<double> magnetization_z;
} // namespace fdtd::materials

fdtd::materials::ModelSpec fdtd::materials::the_model_spec;
fdtd::materials::DispersiveState fdtd::materials::the_dispersive_state;
fdtd::materials::ModelParameters fdtd::materials::the_model_parameters;
fdtd::array1<int> fdtd::materials::background_y;
