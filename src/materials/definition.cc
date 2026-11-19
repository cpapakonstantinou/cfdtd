/**
 * @file materials/definition.cc
 * @brief A material as data, materials/definition.h
 */
#include "materials/definition.h"

#include <charconv>
#include <map>
#include <set>
#include <stdexcept>

namespace fdtd::materials
{
namespace
{
constexpr unsigned bit(Kind k) { return 1u << unsigned(k); }

constexpr unsigned ALL = (1u << kind_count) - 1;
constexpr unsigned DISPERSIVE = ALL & ~bit(Kind::dielectric);
constexpr unsigned FREE_ELECTRONS = bit(Kind::drude_cp) | bit(Kind::drude) | bit(Kind::hydrodynamic) | bit(Kind::drude_lorentz);
constexpr unsigned LORENTZ = bit(Kind::lorentz) | bit(Kind::magnetic_lorentz);
constexpr unsigned CRITICAL_POINTS = bit(Kind::drude_cp) | bit(Kind::hydrodynamic);
constexpr unsigned POLES = CRITICAL_POINTS | bit(Kind::modified_lorentz) | bit(Kind::drude_lorentz);

/// a number of the material section: its key, where it goes, the models that take it
struct Key
{
	std::string_view name;
	double Definition::*member;
	unsigned kinds;
};

constexpr Key KEYS[] = {
	{"n", &Definition::n, bit(Kind::dielectric)},
	{"eps", &Definition::eps, bit(Kind::dielectric)},
	{"conductivity", &Definition::conductivity, ALL},
	{"eps_inf", &Definition::eps_inf, DISPERSIVE},
	{"eps_static", &Definition::eps_static, LORENTZ},
	{"plasma_frequency", &Definition::plasma_frequency, FREE_ELECTRONS},
	{"plasma_damping", &Definition::plasma_damping, FREE_ELECTRONS},
	{"resonance_frequency", &Definition::resonance_frequency, LORENTZ},
	{"damping", &Definition::damping, LORENTZ},
	{"magnetic_coupling", &Definition::magnetic_coupling, bit(Kind::magnetic_lorentz)},
	{"carrier_density", &Definition::carrier_density, bit(Kind::drude)},
	{"effective_mass", &Definition::effective_mass, bit(Kind::drude) | bit(Kind::hydrodynamic)},
	{"electron_density", &Definition::electron_density, FREE_ELECTRONS},
	{"fermi_velocity", &Definition::fermi_velocity, FREE_ELECTRONS},
	{"diffusion", &Definition::diffusion, FREE_ELECTRONS},
	{"sheet_thickness", &Definition::sheet_thickness, bit(Kind::drude_cp)},
	{"sheet_scale", &Definition::sheet_scale, bit(Kind::drude_cp)},
};

/// a number of a pole section
struct PoleKey
{
	std::string_view name;
	double Pole::*member;
	unsigned kinds;
};

constexpr PoleKey POLE_KEYS[] = {
	{"amplitude", &Pole::amplitude, POLES},
	{"frequency", &Pole::frequency, POLES},
	{"damping", &Pole::damping, POLES},
	{"phase", &Pole::phase, CRITICAL_POINTS},
	{"numerator", &Pole::numerator, bit(Kind::modified_lorentz)},
};

/// the keys of the material section that are not numbers
const std::set<std::string> TEXT_KEYS = {"model", "extends", "reference", "chi2", "chi3"};

const std::vector<std::string> COMPONENTS_2 = {
	"d11", "d12", "d13", "d14", "d15", "d16",
	"d21", "d22", "d23", "d24", "d25", "d26",
	"d31", "d32", "d33", "d34", "d35", "d36"};

const std::vector<std::string> COMPONENTS_3 = {
	"xxxx", "xxyy", "xyyx", "xyxy", "xxzz", "xzzx", "xzxz",
	"yyyy", "yyxx", "yxxy", "yxyx", "yyzz", "yzzy", "yzyz",
	"zzzz", "zzxx", "zxxz", "zxzx", "zzyy", "zyyz", "zyzy"};

using Values = std::map<std::string, std::string>;

/// the keys of a material and its poles after extends, the last word of each key wins
struct Raw
{
	Values main;
	std::map<int, Values> poles;
	std::map<std::string, std::string> origin; ///< key -> the section it came from, for the errors
};

std::string format(double v)
{
	char buf[64];
	const auto res = std::to_chars(buf, buf + sizeof(buf), v);
	return std::string(buf, res.ptr);
}

bool valid_name(const std::string &name)
{
	if (name.empty())
		return false;
	for (char ch : name)
		if (!((ch >= 'a' && ch <= 'z') || (ch >= '0' && ch <= '9') || ch == '_' || ch == '-'))
			return false;
	return true;
}

void collect(const io::Ini &ini, const std::string &name, Raw &raw, std::vector<std::string> &chain)
{
	if (!valid_name(name))
		throw std::invalid_argument("material name '" + name + "': use lower-case letters, digits, '_' and '-'");
	for (const auto &c : chain)
		if (c == name)
			throw std::invalid_argument("material '" + chain.front() + "': extends goes round in a circle through '" + name + "'");
	chain.push_back(name);
	const std::string sect = material_section(name);
	if (!ini.is_section(sect))
		throw std::invalid_argument("material '" + name + "' is not defined: no [" + sect + "] section in the run file or in the catalog");
	const std::string parent = io::Ini::lower(ini.get(sect, "extends"));
	if (!parent.empty())
		collect(ini, parent, raw, chain);
	for (const auto &key : ini.keys(sect))
		if (key != "extends")
		{
			raw.main[key] = ini.get(sect, key);
			raw.origin[key] = sect;
		}
	// poles: [material.<name>.pole.<i>], and no other section of the material
	const std::string prefix = sect + ".";
	for (const auto &s : ini.sections())
	{
		if (s.rfind(prefix, 0) != 0)
			continue;
		const std::string rest = s.substr(prefix.size());
		int p = 0;
		if (rest.rfind("pole.", 0) == 0)
		{
			const std::string number = rest.substr(5);
			const auto res = std::from_chars(number.data(), number.data() + number.size(), p);
			if (res.ec != std::errc() || res.ptr != number.data() + number.size())
				p = 0;
		}
		if (p < 1 || p > max_poles)
			throw std::invalid_argument("[" + s + "]: a material has [" + prefix + "pole.1] to [" + prefix + "pole." + std::to_string(max_poles) + "] only");
		for (const auto &key : ini.keys(s))
		{
			raw.poles[p][key] = ini.get(s, key);
			raw.origin["pole." + std::to_string(p) + "." + key] = s;
		}
	}
}

Kind kind_of(const std::string &model, const std::string &sect)
{
	for (int k = 0; k < kind_count; ++k)
		if (KIND_NAMES[k] == model)
			return Kind(k);
	std::string all;
	for (const auto &n : KIND_NAMES)
		all += (all.empty() ? "" : ", ") + std::string(n);
	throw cfg::syntax_error(sect, "model", "'" + model + "' is not a model: " + all);
}

void read_tensor(const io::Ini &ini, int order, const std::string &name, double *values)
{
	const std::string sect = tensor_section(order, name);
	if (!ini.is_section(sect))
		throw std::invalid_argument("chi" + std::to_string(order) + " tensor '" + name + "' is not defined: no [" + sect + "] section in the run file or in the catalog");
	for (const auto &key : ini.keys(sect))
	{
		const int i = tensor_component(order, key);
		if (i < 0)
			throw cfg::syntax_error(sect, key, order == 2 ? "chi2 components are d11 ... d36, d_mn in m/V"
														  : "chi3 components are xxxx, xxyy, ... zyzy, in m^2/V^2");
		values[i] = io::Ini::to_double(sect, key, ini.get(sect, key));
	}
}
} // namespace

Definition vacuum()
{
	Definition d;
	d.name = "vacuum";
	d.kind = Kind::dielectric;
	d.eps = 1;
	return d;
}

std::string material_section(const std::string &name) { return "material." + name; }
std::string pole_section(const std::string &name, int pole) { return "material." + name + ".pole." + std::to_string(pole); }
std::string tensor_section(int order, const std::string &name) { return "chi" + std::to_string(order) + "." + name; }

const std::vector<std::string> &tensor_components(int order) { return order == 2 ? COMPONENTS_2 : COMPONENTS_3; }

int tensor_component(int order, const std::string &key)
{
	const auto &names = tensor_components(order);
	for (size_t i = 0; i < names.size(); ++i)
		if (names[i] == key)
			return int(i);
	return -1;
}

std::vector<std::string> references(const io::Ini &ini, const std::string &name)
{
	const std::string sect = material_section(name);
	std::vector<std::string> out;
	if (!ini.is_section(sect))
		return out;
	if (const std::string parent = io::Ini::lower(ini.get(sect, "extends")); !parent.empty())
		out.push_back(material_section(parent));
	for (int order : {2, 3})
		if (const std::string t = io::Ini::lower(ini.get(sect, "chi" + std::to_string(order))); !t.empty())
			out.push_back(tensor_section(order, t));
	return out;
}

Definition definition(const io::Ini &ini, const std::string &name)
{
	if (name == "vacuum")
	{
		if (ini.is_section(material_section(name)))
			throw std::invalid_argument("[material.vacuum]: vacuum is built in and cannot be redefined");
		return vacuum();
	}
	Raw raw;
	std::vector<std::string> chain;
	collect(ini, name, raw, chain);
	const std::string sect = material_section(name);
	Definition d;
	d.name = name;
	if (!raw.main.count("model"))
		throw cfg::syntax_error(sect, "model", "a material needs a model");
	d.kind = kind_of(io::Ini::lower(raw.main["model"]), raw.origin["model"]);
	const unsigned mine = bit(d.kind);
	for (const auto &[key, value] : raw.main)
	{
		if (TEXT_KEYS.count(key))
			continue;
		const Key *k = nullptr;
		for (const auto &c : KEYS)
			if (c.name == key)
				k = &c;
		if (!k || !(k->kinds & mine))
			throw cfg::syntax_error(raw.origin[key], key, "not a key of model " + std::string(KIND_NAMES[int(d.kind)]));
		d.*(k->member) = io::Ini::to_double(raw.origin[key], key, value);
	}
	if (raw.main.count("reference"))
		d.reference = raw.main["reference"];
	for (const auto &[p, values] : raw.poles)
	{
		if (!(POLES & mine))
			throw std::invalid_argument("material '" + name + "': model " + std::string(KIND_NAMES[int(d.kind)]) + " has no poles");
		if (p != int(d.poles.size()) + 1)
			throw std::invalid_argument("material '" + name + "': the poles are numbered 1, 2, ... without a gap");
		Pole pole;
		for (const auto &[key, value] : values)
		{
			const std::string where = raw.origin["pole." + std::to_string(p) + "." + key];
			const PoleKey *k = nullptr;
			for (const auto &c : POLE_KEYS)
				if (c.name == key)
					k = &c;
			if (!k || !(k->kinds & mine))
				throw cfg::syntax_error(where, key, "not a pole key of model " + std::string(KIND_NAMES[int(d.kind)]));
			pole.*(k->member) = io::Ini::to_double(where, key, value);
		}
		d.poles.push_back(pole);
	}
	if (d.kind == Kind::dielectric && d.n == 0 && d.eps == 0)
		throw std::invalid_argument("material '" + name + "': a dielectric needs n or eps");
	if (d.kind != Kind::dielectric && d.eps_inf == 0)
		throw std::invalid_argument("material '" + name + "': model " + std::string(KIND_NAMES[int(d.kind)]) + " needs eps_inf");
	if (d.kind == Kind::drude_lorentz && d.poles.size() > 1)
		throw std::invalid_argument("material '" + name + "': the Drude-Lorentz model has one pole");
	d.chi2 = io::Ini::lower(raw.main.count("chi2") ? raw.main["chi2"] : "");
	d.chi3 = io::Ini::lower(raw.main.count("chi3") ? raw.main["chi3"] : "");
	if (!d.chi2.empty())
		read_tensor(ini, 2, d.chi2, d.chi2_values.data());
	if (!d.chi3.empty())
		read_tensor(ini, 3, d.chi3, d.chi3_values.data());
	return d;
}

std::string to_ini(const Definition &d)
{
	const Definition none;
	const unsigned mine = bit(d.kind);
	std::string out = "\n[" + material_section(d.name) + "]\nmodel = " + std::string(KIND_NAMES[int(d.kind)]) + "\n";
	if (!d.reference.empty())
		out += "reference = " + d.reference + "\n";
	for (const auto &k : KEYS)
		if ((k.kinds & mine) && d.*(k.member) != none.*(k.member))
			out += std::string(k.name) + " = " + format(d.*(k.member)) + "\n";
	if (!d.chi2.empty())
		out += "chi2 = " + d.chi2 + "\n";
	if (!d.chi3.empty())
		out += "chi3 = " + d.chi3 + "\n";
	for (size_t p = 0; p < d.poles.size(); ++p)
	{
		out += "\n[" + pole_section(d.name, int(p) + 1) + "]\n";
		for (const auto &k : POLE_KEYS)
			if (k.kinds & mine)
				out += std::string(k.name) + " = " + format(d.poles[p].*(k.member)) + "\n";
	}
	return out;
}

std::string tensor_to_ini(int order, const std::string &name, const double *values)
{
	std::string out = "\n[" + tensor_section(order, name) + "]\n";
	const auto &names = tensor_components(order);
	for (size_t i = 0; i < names.size(); ++i)
		if (values[i] != 0)
			out += names[i] + " = " + format(values[i]) + "\n";
	return out;
}
} // namespace fdtd::materials
