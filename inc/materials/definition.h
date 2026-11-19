/**
 * @file materials/definition.h
 * @brief A material as data: its model and the model's parameters, read from INI sections
 *
 * A material is a set of sections named after it, in a run file or in a catalog file, docs/Materials.md:
 * @verbatim
      [material.silver]           model, parameters, extends, chi2, chi3, reference
      [material.silver.pole.1]    one term of the model, numbered from 1
      [chi3.au_isotropic]         a nonlinear tensor named by chi2 = / chi3 =
   @endverbatim
 * `extends = <name>` starts from another material and overrides key by key, pole by pole. Each model
 * accepts the keys of the key table, kinds(); any other key is an error.
 */
#ifndef FDTD_MATERIALS_DEFINITION_H
#define FDTD_MATERIALS_DEFINITION_H

#include <array>
#include <string>
#include <string_view>
#include <vector>

#include "io/ini.h"

namespace fdtd::materials
{
/// the update models, in the order of the former code ranges
enum class Kind : unsigned char
{
	dielectric,
	drude_cp,
	drude,
	lorentz,
	modified_lorentz,
	hydrodynamic,
	drude_lorentz,
	magnetic_lorentz,
};

inline constexpr int kind_count = 8;
inline constexpr std::array<std::string_view, kind_count> KIND_NAMES = {
	"dielectric", "drude_cp", "drude", "lorentz", "modified_lorentz", "hydrodynamic", "drude_lorentz", "magnetic_lorentz"};

inline constexpr int max_poles = 4; ///< terms of one material, at most

/// one term of a model: a critical point, a Lorentz pole of the modified Lorentz or Drude-Lorentz model
struct Pole
{
	double amplitude = 0;
	double frequency = 0; ///< [rad/s]
	double damping = 0;   ///< [rad/s]
	double phase = 0;     ///< [rad], critical points
	double numerator = 0; ///< [rad/s], modified Lorentz: the i omega term of the numerator
	bool operator==(const Pole &) const = default;
};

struct Definition
{
	std::string name;
	Kind kind = Kind::dielectric;
	std::string reference; ///< where the values come from, not read by the solver
	double n = 0;            ///< dielectric: refractive index, eps = n^2
	double eps = 0;          ///< dielectric: relative permittivity, when n is not given
	double conductivity = 0; ///< [S/m]
	double eps_inf = 0;
	double eps_static = 0;          ///< Lorentz, magnetic Lorentz
	double plasma_frequency = 0;    ///< [rad/s]
	double plasma_damping = 0;      ///< [rad/s]
	double resonance_frequency = 0; ///< [rad/s], Lorentz, magnetic Lorentz
	double damping = 0;             ///< [rad/s], Lorentz, magnetic Lorentz
	double magnetic_coupling = 0;   ///< magnetic Lorentz: the a3 term
	double carrier_density = 0;     ///< [m^-3] Drude: plasma frequency from the carriers when given
	double effective_mass = 0;      ///< electron masses; Drude carriers, hydrodynamic; 0 = from the density and plasma frequency
	double electron_density = 0;    ///< [m^-3] free electrons, nonlocal and hydrodynamic models
	double fermi_velocity = 0;      ///< [m/s] nonlocal pressure, times [electrons] fermi_velocity
	double diffusion = 0;           ///< [m^2/s] GNOR diffusion
	double sheet_thickness = 0;     ///< [m] a sheet drawn one cell thick: plasma frequency and amplitudes scaled by thickness / dy
	double sheet_scale = 1;         ///< factor on that scaling
	std::vector<Pole> poles;
	std::string chi2; ///< tensor names, [chi2.<name>] / [chi3.<name>]; empty = linear
	std::string chi3; ///< tensor names, [chi2.<name>] / [chi3.<name>]; empty = linear
	std::array<double, 18> chi2_values{}; ///< d11 .. d36 [m/V] of chi2, as given
	std::array<double, 21> chi3_values{}; ///< xxxx .. zyzy [m^2/V^2] of chi3
	bool operator==(const Definition &) const = default;
};

/// the material every run has, slot 0: the background of every cell nothing is drawn in
Definition vacuum();

/// section names of material `name`: the material, a pole, a tensor
std::string material_section(const std::string &name);
std::string pole_section(const std::string &name, int pole);
std::string tensor_section(int order, const std::string &name);

/// component names of a tensor: d11 .. d36, order 2, or xxxx .. zyzy, order 3
const std::vector<std::string> &tensor_components(int order);
/// position of `key` in tensor_components(order), or -1
int tensor_component(int order, const std::string &key);

/// the sections material `name` needs beyond its own: the material it extends, its tensors
std::vector<std::string> references(const io::Ini &ini, const std::string &name);

/// material `name` from the INI, extends and tensors resolved; vacuum when the INI does not define
/// it and the name is vacuum. \throws syntax_error, std::invalid_argument
Definition definition(const io::Ini &ini, const std::string &name);

/// the sections of a definition, extends resolved: the material and its poles, without the tensors
std::string to_ini(const Definition &d);
/// the section of tensor `name` of order 2 or 3
std::string tensor_to_ini(int order, const std::string &name, const double *values);
} // namespace fdtd::materials

#endif
