/** \file 	config/ini.cc
 * \brief 	implementation for configuration utilties
 * \author 	c. papakonstantinou
 * \date 	March 2023 
 **/
#include <fstream>
#include <errno.h>
#include <cctype>
#include <sstream>
#include <cstdlib>
#include <ctype.h>

#include "config/ini.h"

using namespace fdtd::config;

Config::Config(const std::string& filename)
{
	load(filename);
}

static
void
strip_comment(std::string& line)
{
	std::string::size_type idx = line.find_first_of(";#");

	if (idx != std::string::npos)
		line.erase(idx);
}

/// Section names and keys are case-insensitive: normalize them to lower case.
static
std::string
lower(std::string s)
{
	for (auto& c : s)
		c = static_cast<char>(tolower(static_cast<unsigned char>(c)));
	return s;
}

bool
Config::add_section(const std::string& name)
{
	const bool added = _sections.insert(make_pair(lower(name), keys_type{})).second;
	if (added)
		_order.push_back(lower(name));
	return added;
}

std::vector<std::string> 
Config::sections(void) const 
{
	std::vector<std::string> ret;
	for (const auto& p : _sections)
		ret.push_back(p.first);
	return ret;
}

std::vector<std::string>
Config::sections_in_file_order(void) const
{
	return _order;
}

key_list
Config::keys(const std::string& sect) const
{
	auto cur = _sections.find(lower(sect));
	if (cur == _sections.end())
		return key_list{};

	key_list ret{};
	for (const auto& i : cur->second)
		ret.push_back(i.first);

	return ret;
}

bool 
Config::is_section(const std::string& sect) const 
{
	return _sections.find(lower(sect)) != _sections.end();
}

bool
Config::is_key(const std::string& sect, const std::string& key) const
{
	auto s = _sections.find(lower(sect));
	if (s == _sections.end())
		return false;

	auto k = s->second.find(lower(key));
	if (k == s->second.end() || k->second.empty())
		return false;

	return true;
}

std::string
Config::get(const std::string& sect, const std::string& key, const std::string& dflt) const
{
	auto cur_sect = _sections.find(lower(sect));
	if (cur_sect == _sections.end())
		return dflt;

	auto cur_key = cur_sect->second.find(lower(key));
	if (cur_key == cur_sect->second.end())
		return dflt;

	if (cur_key->second.empty())
		return dflt;

	return cur_key->second.back();
}

value_list
Config::get_all(const std::string& sect, const std::string& key) const 
{
	value_list ret;

	auto cur_sect = _sections.find(lower(sect));
	if (cur_sect == _sections.end())
		return ret;

	auto cur_key = cur_sect->second.find(lower(key));
	if (cur_key == cur_sect->second.end())
		return ret;

	return cur_key->second;
}

unsigned 
Config::get_unsigned(const std::string& sect, const std::string& key, unsigned dflt, unsigned base) const 
{
	std::string raw = get(sect, key);
	if (raw.empty())
		return dflt;
	return to_unsigned(sect, key, raw, base);
}

double 
Config::get_double(const std::string& sect, const std::string& key, double dflt) const 
{
	std::string raw = get(sect, key);
	if (raw.empty())
		return dflt;
	return to_double(sect, key, raw);
}

int 
Config::get_int(const std::string& sect, const std::string& key, int dflt) const 
{
	std::string raw = get(sect, key);
	if (raw.empty())
		return dflt;
	return to_integer(sect, key, raw, 10);
}

bool 
Config::get_bool(const std::string& sect, const std::string& key, bool dflt) const 
{
	std::string raw = lower(get(sect, key));
	if (raw.empty())
		return dflt;
	if (raw == "true" || raw == "yes" || raw == "on")
		return true;
	if (raw == "false" || raw == "no" || raw == "off")
		return false;
	return to_integer(sect, key, raw, 10) != 0;
}

void
Config::trim(std::string& line)
{
	if (line.empty())
		return;

	// start by trimming trailing whitespace
	std::string::size_type idx;
	for (idx = line.size() - 1; idx > 0; --idx)
		if (!isspace(line[idx]))
			break;

	if (idx < line.size() - 1)
		line.erase(idx + 1);

	// now, trim leading whitespace
	for (idx = 0; idx < line.size(); ++idx)
		if (!isspace(line[idx]))
			break;

	if (idx < line.size()) // any non-space characters?
		line = std::string(line, idx);
	else
		line = std::string{}; // empty string
}

bool
Config::split_key_value(const std::string& line, std::string& key, std::string& val, char delim)
{
	std::string::size_type idx = line.find(delim);
	if (idx == std::string::npos)
		return false;

	key = line.substr(0, idx);
	val = line.substr(idx + 1);

	return true;
}

void
Config::load(const std::string& filename)
{
	std::ifstream strm(filename);

	if (!strm)
		throw file_error(filename, errno);

	parse(strm, filename);
}

void
Config::load_string(const std::string& text, const std::string& name)
{
	std::istringstream strm(text);
	parse(strm, name);
}

void
Config::parse(std::istream& strm, const std::string& filename)
{
	std::string line;
	unsigned line_num = 0;
	std::string cur_section;

	while (getline(strm, line))
	{
		++line_num;

		strip_comment(line);
		trim(line);
		if (line.empty())
			continue;

		if (line.front() == '[') // start of section header?
		{
			if (line.back() != ']')
				throw syntax_error{filename, line_num, "section header not terminated with ']'"};
			
			cur_section = line.substr(1, line.size() - 2);
			trim(cur_section);

			if (!add_section(cur_section))
				throw syntax_error{filename, line_num, "duplicate section header " + cur_section};
		}
		else // key-value pair
		{
			std::string key;
			std::string val;

			if (cur_section.empty())
				throw syntax_error{filename, line_num, "key-value pair outside of a [section]"};

			if (!split_key_value(line, key, val))
				throw syntax_error{filename, line_num};

			trim(key);
			if (key.empty())
				throw syntax_error{filename, line_num, "missing key"};

			trim(val);
			add_key_val(cur_section, key, val);
		}
	}
}

void
Config::add_key_val(const std::string& sect, const std::string& key, const std::string& val)
{
	add_key_val(sect, key, value_list{val});
}

void
Config::add_key_val(const std::string& sect, const std::string& key, const value_list& vals)
{
	const std::string lsect = lower(sect);
	const std::string lkey = lower(key);
	auto cur_sect = _sections.find(lsect);
	if (cur_sect == _sections.end())
	{
		add_section(lsect);
		cur_sect = _sections.find(lsect);
	}

	keys_type& section_keys(cur_sect->second);
	// try to drop in a blank value list, if there's nothing there already
	section_keys.insert(make_pair(lkey, value_list{}));
	value_list& cur_key_vals =
		section_keys.find(lkey)->second;
	for (const auto& v : vals)
		cur_key_vals.push_back(v);
}

void
Config::delete_key(const std::string& sect, const std::string& key)
{
	auto cur_sect = _sections.find(lower(sect));
	if (cur_sect == _sections.end())
		return; // nothing to do

	auto cur_key = cur_sect->second.find(lower(key));
	if (cur_key == cur_sect->second.end())
		return; // nothing to do
	cur_key->second = value_list{};
}

void
Config::override(const std::string& sect, const override_list& overrides)
{
	for (auto ov : overrides)
	{
		trim(ov.first);
		if (ov.first.empty())
			continue;
		trim(ov.second);
		add_key_val(sect, ov.first, ov.second);
	}
}

std::string
Config::to_string(void) const
{
	std::ostringstream ss;

	for (auto& sect : _sections)
	{
		ss << "[" << sect.first << "]\n";
		for (auto& key : sect.second)
			for (auto& val : key.second)
				ss << key.first << " = " << val << "\n";
	}

	return ss.str();
}

unsigned
Config::to_unsigned(const std::string& sect, const std::string& key, const std::string& raw, unsigned base)
{
	char* eptr;

	unsigned ret = strtoul(raw.c_str(), &eptr, base);

	if (eptr != raw.c_str() && (*eptr == 0 || isspace(*eptr)))
		return ret;

	throw syntax_error{sect, key, "invalid value " + raw};
}

int
Config::to_integer(const std::string& sect, const std::string& key, const std::string& raw, unsigned base)
{
	char* eptr;

	long ret = strtol(raw.c_str(), &eptr, base);

	if (eptr != raw.c_str() && (*eptr == 0 || isspace(*eptr)))
		return static_cast<int>(ret);

	throw syntax_error{sect, key, "invalid integer value " + raw};
}

double
Config::to_double(const std::string& sect, const std::string& key, const std::string& raw)
{
	char* eptr;

	double ret = strtod(raw.c_str(), &eptr);
	if (eptr != raw.c_str() && (*eptr == 0 || isspace(*eptr)))
		return ret;

	throw syntax_error{sect, key, "invalid value " + raw};
}

std::string
Config::unquote(const std::string& str)
{
	std::string ret;

	ret.reserve(str.size());
	for (std::string::size_type i = 0; i < str.size(); ++i)
		if (str[i] != '"')
			ret.push_back(str[i]);

	return ret;
}
