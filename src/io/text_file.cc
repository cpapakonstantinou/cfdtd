/** \file io/text_file.cc
 *  \brief Line reader shared by the remaining text-file inputs, geometry, temperature, mode files. */
#include "io/text.h"
#include <fstream>
#include <string>

/**
 * @brief Parse a line from a text file
 *
 * @param pf address of fstream
 *
 * @return string which preceeds the whitespace (" ") character
 *
 * @details Reads just a single line, breaks the line by a whitespace, and
 * returns the first portion as a string.
 */
std::string fdtd::io::read_text_entry(std::fstream& pf)
{
	std::string linestring, strd;
	size_t npos;

	getline(pf, linestring);
	// an empty space, translated from Italian
	npos = linestring.find(" ");
	strd = linestring.substr(0, npos);
	return strd;
}
