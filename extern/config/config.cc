/** \file 	config.cc
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

#include "config.h"

using namespace std;
using namespace cfg;


Config::Config(const string& filename)
{
	load(filename);
}

static
void
strip_comment(string& line)
{
	string::size_type idx = line.find(';');

	if (idx != string::npos)
		line.erase(idx);
}

bool 
Config::add_section(const string& name) 
{
	return sections_.insert(make_pair(name, keys_type{})).second;
}

vector<string> 
Config::sections(void) const 
{
	vector<string> ret;
	for (const auto& p : sections_)
		ret.push_back(p.first);
	return ret;
}

key_list 
Config::keys(const string& sect) const
{
	auto cur = sections_.find(sect);
	if (cur == sections_.end())
		return key_list{};

	key_list ret{};
	for (const auto& i : cur->second)
		ret.push_back(i.first);

	return ret;
}

bool 
Config::is_section(const string& sect) const 
{
	return sections_.find(sect) != sections_.end();
}

bool 
Config::is_key(const string& sect, const string& key) const 
{
	auto s = sections_.find(sect);
	if (s == sections_.end())
		return false;

	if (s->second.find(key) == s->second.end())
		return false;

	return true;
}

string
Config::get(const string& sect, const string& key, const string& dflt) const
{
	auto cur_sect = sections_.find(sect);
	if (cur_sect == sections_.end())
		return dflt;

	auto cur_key = cur_sect->second.find(key);
	if (cur_key == cur_sect->second.end())
		return dflt;

	if (cur_key->second.empty())
		return dflt;

	return cur_key->second.back();
}

value_list
Config::get_all(const string& sect, const string& key) const 
{
	value_list ret;

	auto cur_sect = sections_.find(sect);
	if (cur_sect == sections_.end())
		return ret;

	auto cur_key = cur_sect->second.find(key);
	if (cur_key == cur_sect->second.end())
		return ret;

	return cur_key->second;
}

unsigned 
Config::get_unsigned(const string& sect, const string& key, unsigned dflt, unsigned base) const 
{
	string raw = get(sect, key);
	if (raw.empty())
		return dflt;
	return to_unsigned(sect, key, raw, base);
}


double 
Config::get_double(const string& sect, const string& key, double dflt) const 
{
	string raw = get(sect, key);
	if (raw.empty())
		return dflt;
	return to_double(sect, key, raw);
}

int 
Config::get_int(const string& sect, const string& key, int dflt) const 
{
	string raw = get(sect, key);
	if (raw.empty())
		return dflt;
	return to_integer(sect, key, raw, 10);
}


bool 
Config::get_bool(const string& sect, const string& key, bool dflt) const 
{
	string raw = get(sect, key);
	if (raw.empty())
		return dflt;
	int c = to_integer(sect, key, raw, 10);
	if (c)
		return true;
	else
		return false;
}


void
Config::trim(string& line)
{
	if (line.empty())
		return;

	// start by trimming trailing whitespace
	string::size_type idx;
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
		line = string(line, idx);
	else
		line = string{}; // empty string
}

bool
Config::split_key_value(const string& line, string& key, string& val, char delim)
{
	string::size_type idx = line.find(delim);
	if (idx == string::npos)
		return false;

	key = line.substr(0, idx);
	val = line.substr(idx + 1);

	return true;
}


void
Config::load(const string& filename)
{
	ifstream strm(filename);

	if (!strm)
		throw file_error(filename, errno);

	string line;
	unsigned line_num = 0;
	string cur_section;

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
			
			if (!add_section(cur_section))
				throw syntax_error{filename, line_num, "duplicate section header " + cur_section};
		}
		else // key-value pair
		{
			string key, val;

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
Config::add_key_val(const string& sect, const string& key, const string& val)
{
	add_key_val(sect, key, value_list{val});
}

void
Config::add_key_val(const string& sect, const string& key, const value_list& vals)
{
	auto cur_sect = sections_.find(sect);
	if (cur_sect == sections_.end())
	{
		add_section(sect);
		cur_sect = sections_.find(sect);
	}

	keys_type& section_keys(cur_sect->second);
	// try to drop in a blank value list, if there's nothing there already
	section_keys.insert(make_pair(key, value_list{}));
	value_list& cur_key_vals = 
		section_keys.find(key)->second;
	for (const auto& v : vals)
		cur_key_vals.push_back(v);
}


void
Config::delete_key(const string& sect, const string& key)
{
	auto cur_sect = sections_.find(sect);
	if (cur_sect == sections_.end())
		return; // nothing to do

	auto cur_key = cur_sect->second.find(key);
	if (cur_key == cur_sect->second.end())
		return; // nothing to do
	cur_key->second = value_list{};
}


void
Config::override(const string& sect, const override_list& overrides)
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

string
Config::to_string(void) const
{
	ostringstream ss;

	for (auto& sect : sections_)
	{
		ss << "[" << sect.first << "]\n";
		for (auto& key : sect.second)
			for (auto& val : key.second)
				ss << key.first << " = " << val << "\n";
	}

	return ss.str();
}


unsigned
Config::to_unsigned(const string& sect, const string& key, const string& raw, unsigned base)
{
	char* eptr;

	unsigned ret = strtoul(raw.c_str(), &eptr, base);

	if (eptr != raw.c_str() && (*eptr == 0 || isspace(*eptr)))
		return ret;

	throw syntax_error{sect, key, "invalid value " + raw};
}

int 
Config::to_integer(const string& sect, const string& key, const string& raw, unsigned base)
{
	return static_cast<int>(to_unsigned(sect, key, raw, base));
}

double
Config::to_double(const string& sect, const string& key, const string& raw)
{
	char* eptr;

	double ret = strtod(raw.c_str(), &eptr);
	if (eptr != raw.c_str() && (*eptr == 0 || isspace(*eptr)))
		return ret;

	throw syntax_error{sect, key, "invalid value " + raw};
}

string
Config::unquote(const string& str)
{
	string ret;

	ret.reserve(str.size());
	for (string::size_type i = 0; i < str.size(); ++i)
		if (str[i] != '"')
			ret.push_back(str[i]);

	return ret;
}
