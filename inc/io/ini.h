/** \file   io/ini.h
 *  \brief  The solver's INI dialect on cfg::Config, extern/config.
 *
 *  cfg::Config stores and queries the sections; Ini reads them, from text so that the
 *  configuration can be broadcast over MPI, with these rules:
 *  - section names and keys are case-insensitive, values keep their case;
 *  - `;` and `#` start a comment;
 *  - a key outside a section, or a section named twice, is an error;
 *  - booleans are true/yes/on/false/no/off or an integer; integers may be negative;
 *  - the sections are also listed in the order they first appear.
 *  What cfg::Config would need to offer instead is listed in extern/README.md.
 */
#ifndef FDTD_IO_INI_H
#define FDTD_IO_INI_H

#include <string>
#include <vector>

#include <config.h>

namespace fdtd::io
{
class Ini : public cfg::Config
{
public:
	/// parse INI text; `name` labels the errors
	void load_string(const std::string& text, const std::string& name = "<memory>");

	/// add values to a section, as cfg::Config::override, with case-insensitive names
	void override(const std::string& sect, const override_list& overrides);

	bool get_bool(const std::string& sect, const std::string& key, bool dflt = false) const;
	int get_int(const std::string& sect, const std::string& key, int dflt = 0) const;

	/// a signed integer value
	static int to_integer(const std::string& sect, const std::string& key, const std::string& raw);

	/// section names in the order they first appeared
	const std::vector<std::string>& sections_in_file_order() const { return _order; }

	static std::string lower(std::string s);

private:
	void add_section_name(const std::string& name);

	std::vector<std::string> _order;
};
} // namespace fdtd::io

#endif // FDTD_IO_INI_H
