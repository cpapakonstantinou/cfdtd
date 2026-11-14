/**
 * @file app/archive.h
 * @brief A binary archive that saves or loads through the same calls
 *
 * A module describes its state once (`ar.array(E); ar.array(H);`); the archive writes it when saving
 * and reads it back when loading, so the two can not diverge. Arrays are loaded into arrays that are
 * already allocated: the stored element count is checked against the file.
 */
#ifndef FDTD_APP_ARCHIVE_H
#define FDTD_APP_ARCHIVE_H

#include "core/carray.h"
#include "materials/material_field.h"
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace fdtd::io
{
class Archive
{
public:
	Archive(const std::string &path, bool saving) : _path(path), _saving(saving), _file(std::fopen(path.c_str(), saving ? "wb" : "rb"))
	{
		if (!_file)
			throw std::runtime_error("cannot open state file " + path);
	}
	~Archive()
	{
		if (_file)
			std::fclose(_file);
	}
	Archive(const Archive &) = delete;
	Archive &operator=(const Archive &) = delete;

	bool saving() const { return _saving; }

	/// a value that is part of the state
	template <class T>
	void value(T &v)
	{
		static_assert(std::is_trivially_copyable_v<T>);
		raw(&v, sizeof v);
	}
	/// a value that must be the same in the file and in this run
	template <class T>
	void same(T v, const char *what)
	{
		T in_file = v;
		value(in_file);
		if (!_saving && in_file != v)
			throw std::runtime_error("state file " + _path + " does not match this run: " + what);
	}
	/// the contents of an allocated array, nothing for an empty one
	template <class T, std::size_t N, std::size_t A>
	void array(carray<T, N, A> &a)
	{
		if (a && a.size() > 0 && a.stored() == 0)
			throw std::logic_error("state archive: an array without a stored count");
		same(std::uint64_t(a.stored()), "array size");
		if (a.stored())
			raw(a.storage(), a.stored() * sizeof(T));
	}
	/// every element of a C array of arrays
	template <class T, std::size_t K>
	void array(T (&a)[K])
	{
		for (T &x : a)
			array(x);
	}
	/// the fields of every material present
	template <class T, std::size_t M>
	void array(material_field<T, M> &f)
	{
		f.for_each([&](ctensor<T> &a) { array(a); });
	}

private:
	void raw(void *data, std::size_t bytes)
	{
		const std::size_t done = _saving ? std::fwrite(data, 1, bytes, _file) : std::fread(data, 1, bytes, _file);
		if (done != bytes)
			throw std::runtime_error(std::string(_saving ? "cannot write state file " : "state file is too short: ") + _path);
	}
	std::string _path;
	bool _saving;
	std::FILE *_file;
};
} // namespace fdtd::io

#endif
