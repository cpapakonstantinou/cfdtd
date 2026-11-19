/** \file   io/ini.cc
 *  \brief  The solver's INI dialect on cfg::Config.
 */
#include "io/ini.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <sstream>

#include <exceptions.h>

using fdtd::io::Ini;

std::string Ini::lower(std::string s)
{
	std::ranges::transform(s, s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	return s;
}

void Ini::add_section_name(const std::string& name)
{
	add_section(name);
	_order.push_back(name);
}

void Ini::load_string(const std::string& text, const std::string& name)
{
	std::istringstream strm(text);
	std::string line;
	unsigned line_num = 0;
	std::string section;
	while (std::getline(strm, line))
	{
		++line_num;
		line = line.substr(0, line.find_first_of(";#"));
		trim(line);
		if (line.empty())
			continue;
		if (line.front() == '[')
		{
			if (line.back() != ']')
				throw cfg::syntax_error{name, line_num, "section header not terminated with ']'"};
			section = line.substr(1, line.size() - 2);
			trim(section);
			section = lower(section);
			if (is_section(section))
				throw cfg::syntax_error{name, line_num, "duplicate section header " + section};
			add_section_name(section);
			continue;
		}
		std::string key;
		std::string val;
		if (!split_key_value(line, key, val))
			throw cfg::syntax_error{name, line_num};
		trim(key);
		if (key.empty())
			throw cfg::syntax_error{name, line_num, "missing key"};
		if (section.empty())
			throw cfg::syntax_error{name, line_num, "key-value pair outside of a [section]"};
		trim(val);
		add_key_val(section, lower(key), val);
	}
}

void Ini::override(const std::string& sect, const override_list& overrides)
{
	const std::string s = lower(sect);
	if (!is_section(s))
		add_section_name(s);
	override_list folded;
	for (const auto& [key, val] : overrides)
		folded.emplace_back(lower(key), val);
	Config::override(s, folded);
}

bool Ini::get_bool(const std::string& sect, const std::string& key, bool dflt) const
{
	const std::string raw = lower(get(sect, key));
	if (raw.empty())
		return dflt;
	if (raw == "true" || raw == "yes" || raw == "on")
		return true;
	if (raw == "false" || raw == "no" || raw == "off")
		return false;
	return to_integer(sect, key, raw) != 0;
}

int Ini::get_int(const std::string& sect, const std::string& key, int dflt) const
{
	const std::string raw = get(sect, key);
	return raw.empty() ? dflt : to_integer(sect, key, raw);
}

int Ini::to_integer(const std::string& sect, const std::string& key, const std::string& raw)
{
	char* end = nullptr;
	const long v = std::strtol(raw.c_str(), &end, 10);
	if (end != raw.c_str() && (*end == 0 || std::isspace(static_cast<unsigned char>(*end))))
		return static_cast<int>(v);
	throw cfg::syntax_error{sect, key, "invalid integer value " + raw};
}
