/** \file   config/chi_tensors.cc
 *  \brief  Component names and built-in presets of the chi2 / chi3 tensors, config/chi_tensors.h.
 *
 *  The presets are the tensors formerly shipped as INIDEF/chi2/ and INIDEF/chi3/ files, with the same values.
 *
 *  References: Boyd, Nonlinear Optics, 3rd ed., 2007; Butcher & Cotter, The Elements of Nonlinear Optics, 1990;
 *  Sutherland, Handbook of Nonlinear Optics, 2003.
 **/
#include "config/chi_tensors.h"

#include <cctype>
#include <utility>

namespace fdtd::config
{
namespace
{

const std::vector<std::string> COMPONENTS_2 = {
	"d11", "d12", "d13", "d14", "d15", "d16",
	"d21", "d22", "d23", "d24", "d25", "d26",
	"d31", "d32", "d33", "d34", "d35", "d36"};

const std::vector<std::string> COMPONENTS_3 = {
	"xxxx", "xxyy", "xyyx", "xyxy", "xxzz", "xzzx", "xzxz",
	"yyyy", "yyxx", "yxxy", "yxyx", "yyzz", "yzzy", "yzyz",
	"zzzz", "zzxx", "zxxz", "zxzx", "zzyy", "zyyz", "zyzy"};

ChiTensor tensor(int order, const char* name, std::initializer_list<std::pair<const char*, double>> components)
{
	ChiTensor t;
	t.name = name;
	t.order = order;
	for (const auto& [key, value] : components)
		t.values[chi_component_index(order, key)] = value;
	return t;
}

/// chi3 of an isotropic or cubic (m3m) medium: xxxx = yyyy = zzzz = `diagonal`, every other component `mixed`.
ChiTensor isotropic(const char* name, double diagonal, double mixed)
{
	ChiTensor t = tensor(3, name, {});
	for (int i = 0; i < 21; ++i)
		t.values[i] = (i % 7 == 0) ? diagonal : mixed;
	return t;
}

const std::vector<ChiTensor> PRESETS_2 = {
	tensor(2, "Zero", {}),                                                           // centrosymmetric, m3m, isotropic media
	tensor(2, "BBO", {{"d11", 1.6e-12}, {"d15", 0.15e-12}, {"d22", 2.1e-12}, {"d31", 0.15e-12}}), // 3m; Sutherland p. 296
	tensor(2, "LiNbO3", {{"d15", -5.95e-12}, {"d31", -5.95e-12}, {"d33", -34.4e-12}}),  // 3m; Sutherland p. 300
	tensor(2, "Diag", {{"d11", 1e-12}, {"d22", 1e-12}, {"d33", 1e-12}}),               // test tensor
};

const std::vector<ChiTensor> PRESETS_3 = {
	tensor(3, "Zero", {}),
	isotropic("Air", 1.7e-25, 0.567e-25),
	isotropic("Au", 7.6e-19, 2.533e-19),     // m3m: xxxx = xxyy + xyxy + xyyx
	tensor(3, "AuDiag", {{"xxxx", 7.6e-19}, {"yyyy", 7.6e-19}, {"zzzz", 7.6e-19}}),
	isotropic("SiO2", 2.5e-22, 0.833e-22),   // fused silica
	tensor(3, "SiO2Diag", {{"xxxx", 2.5e-22}, {"yyyy", 2.5e-22}, {"zzzz", 2.5e-22}}),
};

bool same_name(const std::string& a, const std::string& b)
{
	if (a.size() != b.size())
		return false;
	for (size_t i = 0; i < a.size(); ++i)
		if (std::tolower(static_cast<unsigned char>(a[i])) != std::tolower(static_cast<unsigned char>(b[i])))
			return false;
	return true;
}

} // namespace

const std::vector<std::string>& chi_components(int order) { return order == 2 ? COMPONENTS_2 : COMPONENTS_3; }

int chi_component_index(int order, const std::string& key)
{
	const auto& names = chi_components(order);
	for (size_t i = 0; i < names.size(); ++i)
		if (same_name(names[i], key))
			return int(i);
	return -1;
}

const std::vector<ChiTensor>& chi_presets(int order) { return order == 2 ? PRESETS_2 : PRESETS_3; }

std::vector<std::string> chi_preset_names(int order)
{
	std::vector<std::string> names;
	for (const auto& t : chi_presets(order))
		names.push_back(t.name);
	return names;
}

const ChiTensor* find_chi_tensor(int order, const std::string& name, const std::vector<ChiTensor>& user)
{
	for (const auto& t : user)
		if (t.order == order && same_name(t.name, name))
			return &t;
	for (const auto& t : chi_presets(order))
		if (same_name(t.name, name))
			return &t;
	return nullptr;
}

} // namespace fdtd::config
