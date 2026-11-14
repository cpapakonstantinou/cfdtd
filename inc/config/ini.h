#ifndef __CONFIG_H__
#define __CONFIG_H__

/** \file 	config/ini.h
 * \brief 	definition for configuration utilties
 * \author 	c. papakonstantinou
 * \date 	March 2023 
 **/
#include <string>
#include <cstring>
#include <istream>
#include <list>
#include <map>
#include <vector>
#include "config/exceptions.h"

/** Configuration file parser and utilities.
 *
 * Create an instance of this class, pointing it at a configuration file, to
 * parse it.
 *
 * The configuration file follows standard Windows "INI" file syntax.
 * Comments are prefixed with a semicolon.  Section names are enclosed
 * in square brackets.  Each line in the file that is not a section
 * name is a key=value pair.  Note that all key-value pairs must be
 * enclosed within a section -- it is a syntax error if it is not.
 * Leading and trailing whitespace in both keys and values are
 * stripped.
 *
 * Section names must all be unique.  However, within a section,
 * multiple entries with the same key may be present.
 *
 * Section names and keys are not case sensitive; however, case is
 * preserved for values.
 *
 * Additional external configuration values may be passed in via the
 * override() method.  This is useful for additional configuration
 * options which are passed in via command-line parameters or some
 * other mechanism.
 */
namespace fdtd::config
{
	/// Type of a single value associated with a key.
	using value_type = std::string;
	/// All of the values associated with a key.
	using value_list = std::list<value_type>;
	/// A list of keys.
	using key_list = std::list<std::string>;

	class Config 
	{
		public:

		private:
			/// Mapping between key names and their values
			using keys_type = std::map<std::string, value_list>;
			/// Mapping between section and their keys, and values
			using sections_type = std::map<std::string, keys_type>;

			sections_type _sections; ///< Known sections.
			std::vector<std::string> _order; ///< Section names in the order they first appeared.

		protected:
			/** \brief Add a key-value pair to the specified section.
			*
			* \param sect Section to add to.
			*
			* \param key Key name to add.
			*
			* \param val Value to add.
			*
			*/
			void add_key_val(const std::string& sect, const std::string& key, const std::string& val);
			/** \brief Add a key-value set to the specified section.
			*
			* \param sect Section to add to.
			*
			* \param key Key name to add to.
			*
			* \param vals Values to add.
			*/
			void add_key_val(const std::string& sect, const std::string& key, const value_list& vals);
			/** \brief Delete all key-value data from a specified section.
			*
			* \param sect Section to delete from.
			*
			* \param key Key to delete.
			*/
			void delete_key(const std::string& sect, const std::string& key);

		public:

			/// Default constructor - empty configuration.
			Config() = default;
			virtual ~Config() = default;

			/** \brief Open and parse a configuration file.
			 *
			 * \param filename Name of the file to parse, relative to cwd.
			 *
			 * \throws file_error, syntax_error.
			 */
			Config(const std::string& filename);

			/** \brief Open and parse a configuration file.
			*
			* \param filename Name of the file to parse, relative to cwd.
			*
			* \throws file_error, syntax_error.
			*/
			void load(const std::string& filename);

			/** \brief Parse configuration text held in memory, e.g. broadcast over MPI.
			*
			* \param text INI text.
			*
			* \param name Name used in error messages.
			*
			* \throws syntax_error.
			*/
			void load_string(const std::string& text, const std::string& name = "<memory>");

			/** \brief Parse configuration text from a stream. \sa load */
			void parse(std::istream& strm, const std::string& name);

			/** \brief Add a section to the configuration file.
			 *
			 * \param name Section name.
			 *
			 * \returns true if the section was addedd successfully, false otherwise.
			 */
			bool add_section(const std::string& name); 

			/** List for configuration overrides.  The first member is the name
			  to override; the second is the value to override with. */
			using override_list = std::list<std::pair<std::string, value_type> >;
			/** \brief Add external overrides to the config file.
			*
			* \param sect Section to add overrides to.  Note that the section
			* will be created if it does not already exist.
			*
			* \param overrides Overrides to include.
			*
			* \note Overrides are treated exactly as if they were lines in the
			* configuration file: if the keys exist already, or if multiple
			* values are supplied for the same key, the values are appended in
			* the order supplied.
			*/
			void override(const std::string &, const override_list &);

			/** \brief Return the available section names. */
			std::vector<std::string> sections(void) const; ///< sorted by name
			std::vector<std::string> sections_in_file_order(void) const; ///< in the order they first appeared

			/** \brief Return the available keys within a section. */
			key_list keys(const std::string& sect) const;

			/** \brief Is a particular section available? 
			 *
			 * \param sect Section name.
			 *
			 * \returns true if the section has been defined. */
			bool is_section(const std::string& sect) const;

			/** \brief Is a particular key in a section available? */
			bool is_key(const std::string& sect, const std::string& key) const;

			/** \brief Convert a config file into a string, for output.
			*
			* \returns The file, in a pretty string form.
			*/
			std::string to_string(void) const;

			/** \brief Return a configuration entry.
			*
			* \param sect Section to search.
			*
			* \param key Name of parameter to search for.
			*
			* \param dflt Default value to return if not present in the config file.
			*
			* \returns The configuration value.  If multiple values with the same
			* name are present, returns the last one.
			*/
			std::string get(const std::string& sect, const std::string& key, const std::string& dflt = std::string("")) const;

			/** \brief Return all configuration entries with a certain key.
			*
			* \param sect Section to search.
			*
			* \param key Name of parameter to search for.
			*
			* \returns The configuration entries, or an empty list if none are
			* found in the config file.  The order of entries matches that of the
			* underlying file.
			*/
			value_list get_all(const std::string& sect, const std::string& key) const;

			/** \brief Get a value as an unsigned quantity.
			 *
			 * \param sect Section to retrieve from.
			 *
			 * \param key Key to retrieve.
			 *
			 * \param dflt Default value to return.
			 *
			 * \param base Base to use for conversions.
			 *
			 * \returns The value.
			 *
			 * \throws syntax_error.
			 */
			unsigned get_unsigned(const std::string& sect, const std::string& key, unsigned dflt = 0, unsigned base = 10) const;

			/** \brief Get a value as a floating-point value.
			 *
			 * \param sect Section to retrieve from.
			 *
			 * \param key Key to retrieve.
			 *
			 * \param dflt Default value to return.
			 *
			 * \returns The value.
			 *
			 * \throws syntax_error.
			 */
			double get_double(const std::string& sect, const std::string& key, double dflt = 0) const;

			/** \brief Get a value as a signed quantity.
			 *
			 * \param sect Section to retrieve from.
			 *
			 * \param key Key to retrieve.
			 *
			 * \param dflt Default value to return.
			 *
			 * \returns The value.
			 *
			 * \throws syntax_error.
			 */
			int get_int(const std::string& sect, const std::string& key, int dflt = 0) const;

			/** \brief Get a value as a boolean quantity.
			 *
			 * \param sect Section to retrieve from.
			 *
			 * \param key Key to retrieve.
			 *
			 * \param dflt Default value to return.
			 *
			 * \returns The value.
			 *
			 * \throws syntax_error.
			 */
			bool get_bool(const std::string& sect, const std::string& key, bool dflt = false) const;

			/** \brief Convert a string value to an unsigned integer, with error
			* checking.
			*
			* \param sect Section key belong to.
			*
			* \param key Key name.
			*
			* \param raw Raw string to convert.
			*
			* \param base Base that string in in; if zero, do auto-detection per C rules.
			*
			* \returns The converted string.
			* 
			* \throws syntax_error.
			*/
			static unsigned to_unsigned(const std::string& sect, const std::string& key, const std::string& raw, unsigned base = 10);

			/** \brief Convert a string value to a signed integer. \sa to_unsigned */
			static int to_integer(const std::string& sect, const std::string& key, const std::string& raw, unsigned base = 10);

			/** \brief Convert a string value to a floating point value, with
			* error checking.
			*
			* \param sect Section key belongs to.
			*
			* \param key Key name.
			*
			* \param raw Raw string to convert.
			*
			* \returns The converted value.
			*
			* \throws syntax_error.
			*/
			static double to_double(const std::string& sect, const std::string& key, const std::string& raw);

			/** \brief Remove quotes from a string value.
			*
			* \param str String to unquote.
			*
			* \returns The unquoted string.
			*/
			static std::string unquote(const std::string& str);

			/** \brief Split a line into a key-value pair.
			*
			* \param line Line to split.
			*
			* \param key Key output.
			*
			* \param val Value output.
			*
			* \param delim Delimiter to use for the split.
			*
			* \returns True if the split was successful, false otherwise.
			*/
			static bool split_key_value(const std::string& line, std::string& key, std::string& val, char delim = '=');

			/** \brief Trim leading and trailing whitespace from a line.
			*
			* \param line Line to trim.
			*/
			static void trim(std::string& line);

	};
}//namespace fdtd::config

#endif //__CONFIG_H__