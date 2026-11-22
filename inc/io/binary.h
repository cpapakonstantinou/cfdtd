#ifndef IO_H
#define IO_H

#include <string>
#include <utility>
#include "core/array.h"
#include "io/h5io.h"
#include "mpi.h"

// Saving Binary Files
template <typename T> void save_as_binary_2d(std::string const &, T ***, T ***, T ***, MPI_Comm const &, int const &,
int* const &, int* const &, int* const &, int* const &, int const &, int const &, bool const &, bool const &, bool const&, const h5io::FieldNote *note = nullptr);

/// `note`: a complex field's scale and its attributes, see h5io::FieldNote
template <typename T> void save_as_binary_3d(std::string const &, T ***, T ***, T ***, MPI_Comm const &, int* const &, int* const &,
int* const &, int* const &, bool const &, bool const &, bool const&, const h5io::FieldNote *note = nullptr);

/// carray owners: forward their views
template <typename T, typename... Rest>
void save_as_binary_2d(std::string const &f, fdtd::array3<T> &x, fdtd::array3<T> &y, fdtd::array3<T> &z, Rest &&...rest)
{
	save_as_binary_2d<T>(f, x.get(), y.get(), z.get(), std::forward<Rest>(rest)...);
}
template <typename T, typename... Rest>
void save_as_binary_3d(std::string const &f, fdtd::array3<T> &x, fdtd::array3<T> &y, fdtd::array3<T> &z, Rest &&...rest)
{
	save_as_binary_3d<T>(f, x.get(), y.get(), z.get(), std::forward<Rest>(rest)...);
}

bool ensure_directory(std::string &dirPath);
bool check_file_exists(const std::string &fileName);

#endif  