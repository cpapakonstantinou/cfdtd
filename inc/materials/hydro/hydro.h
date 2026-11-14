/**
 * @file materials/hydro/hydro.h
 * @brief The electron-gas models: the hydrodynamic model, hydrodynamics.cc, the nonlocal update and
 *        the electron density of the local models (electron_gas.cc)
 *
 * update_p{x,y,z}[_2]: the polarisation update of one component, first and second pass;
 * setup_eps_inf_{x,y,z}: the background permittivity of the component's cells. The per-component
 * copies wait for the physics review, roadmap §5.
 */
#ifndef FDTD_MATERIALS_HYDRO_HYDRO_H
#define FDTD_MATERIALS_HYDRO_HYDRO_H

namespace fdtd::materials::hydro
{
/// hydrodynamics.cc
void update_electron_density();
void update_hydrodynamics();
void update_px();
void update_py();
void update_pz();
void update_px_2();
void update_py_2();
void update_pz_2();
void setup_eps_inf_x();
void setup_eps_inf_y();
void setup_eps_inf_z();
/// electron_gas.cc: the nonlocal, GNOR, polarisation update and the electron density of the local models
void update_nonlocal();
void update_electron_density_local();
} // namespace fdtd::materials::hydro

#endif
