#ifndef IO_H
#define IO_H

#include <string>
#include <utility>
#include "core/carray.h"
#include "mpi.h"

// Saving Binary Files
template <typename T> void save_as_binary_2d(std::string const &, T ***, T ***, T ***, MPI_Comm const &, int const &, 
int* const &, int* const &, int* const &, int* const &, int const &, int const &, bool const &, bool const &, bool const&);

template <typename T> void save_as_binary_3d(std::string const &, T ***, T ***, T ***, MPI_Comm const &, int* const &, int* const &, 
int* const &, int* const &, bool const &, bool const &, bool const&);

/// carray owners: forward their views
template <typename T, typename... Rest>
void save_as_binary_2d(std::string const &f, ctensor<T> &x, ctensor<T> &y, ctensor<T> &z, Rest &&...rest)
{
	save_as_binary_2d<T>(f, x.get(), y.get(), z.get(), std::forward<Rest>(rest)...);
}
template <typename T, typename... Rest>
void save_as_binary_3d(std::string const &f, ctensor<T> &x, ctensor<T> &y, ctensor<T> &z, Rest &&...rest)
{
	save_as_binary_3d<T>(f, x.get(), y.get(), z.get(), std::forward<Rest>(rest)...);
}

bool ensure_directory(std::string &dirPath);
bool check_file_exists(const std::string &fileName);

#endif  