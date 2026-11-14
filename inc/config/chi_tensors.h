/** \file   config/chi_tensors.h
 *  \brief  Nonlinear susceptibility tensors for [nonlinear]: named components, built-in presets, [chi2.<name>] / [chi3.<name>] sections.
 *
 *  A tensor is given in the INI by its non-zero components; the components not listed are zero.
 *
 *  - chi2: the 18 contracted coefficients d11 ... d36 in m/V (chi2_ijk = 2 d_mn).
 *  - chi3: the 21 components xxxx, xxyy, ... zyzy in m^2/V^2.
 *
 *  ```ini
 *  [nonlinear]
 *  shg = antenna
 *  chi_structure = my_crystal   ; a [chi2.my_crystal] section, or a preset, BBO, LiNbO3, ...
 *
 *  [chi2.my_crystal]
 *  d22 = 2.1e-12
 *  d31 = 0.15e-12
 *  ```
 *
 *  A section of the same name as a preset replaces that preset.
 **/
#pragma once

#include <array>
#include <string>
#include <vector>

namespace fdtd::config
{

/** One named tensor: `order` is 2 or 3; `values` holds the components in the order of chi_components(order), 18 or 21 used. */
struct ChiTensor
{
	std::string name;
	int order = 2;
	std::array<double, 21> values{};
	bool operator==(const ChiTensor&) const = default;
};

const std::vector<std::string>& chi_components(int order);   ///< d11 ... d36, order 2, or xxxx ... zyzy, order 3
int chi_component_index(int order, const std::string& key);    ///< position of `key` in chi_components(order), or -1
const std::vector<ChiTensor>& chi_presets(int order);        ///< the built-in tensors
std::vector<std::string> chi_preset_names(int order);

/** The tensor called `name`: a section in `user` first, then the presets; nullptr if neither has it. Case-insensitive. */
const ChiTensor* find_chi_tensor(int order, const std::string& name, const std::vector<ChiTensor>& user);

} // namespace fdtd::config
