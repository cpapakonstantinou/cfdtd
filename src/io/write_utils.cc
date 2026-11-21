/**
 * @file io/write_utils.cc
 * @author FDTD-Lab
 * @brief Functions to write out binary files and check output paths
 * @version 0.1
 * @date 2022-04-12
 * 
 * @copyright Copyright 2022
 */
#include <string>
#include "mpi.h"
#include "io/text.h"
#include <cstdio>
#include <cstdarg>
#include <typeinfo>
#include <type_traits>
#include <filesystem>
#include "io/binary.h"
#include "io/h5io.h"
#include <vector>

using namespace std;
namespace fs = std::filesystem;

/**
 * @brief Checks to see if path exists and is directory
 * If doesn't exist, will create the path
 *
 * @param dirPath [in] path string to directory you wish to check
 * @return true if error eccured
 * @return false if there were no issues
 */
bool ensure_directory(std::string &dirPath)
{
	fs::path tmpPath{ dirPath };

	if (fs::exists(tmpPath))
	{
		if (!fs::is_directory(dirPath))
		{
			std::fprintf(stderr, "error: %s exists and is not a directory\n", dirPath.c_str());
			return true;
		}
		else
			return false;
	}
	return !fs::create_directory(dirPath);
}

/**
 * @brief Checks if a file exists
 *
 * @param fileName [in] path and name of file to check
 * @return true if file found
 * @return false if file not found
 */
bool check_file_exists(const std::string &fileName) {
	if (FILE *tmpFile = fopen(fileName.c_str(), "r"))
	{
		fclose(tmpFile);
		return true;
	} else {
		return false;
	}
}

/**
 * @brief gets the time passed between startTime and end time.
 * @param hours variable to store hours passed
 * @param minutes variable to store minutes passed
 * @param seconds variable to store seconds passed
 */
void get_time(clock_t startTime, clock_t endTime, float &hours, float &minutes, float &seconds)
{
	float secselapsed = (endTime - startTime) / CLOCKS_PER_SEC;
	float r1;
	float r2;

	r1 = 3600. * modf(secselapsed / 3600., &hours);
	r2 = modf(r1 / 60., &minutes);
	seconds = r2 * 60;
}

/**
 * @brief Appends a line to <dir>/run.log with the time since start, creating the directory
 * @return true if the directory could not be made
 */
bool fdtd::io::write_log(const std::string &text, const std::string &dir, clock_t start, bool verbose)
{
	clock_t currentTime = clock();
	float hours = 0;
	float minutes = 0;
	float seconds = 0;
	std::string path = dir;
	bool error = ensure_directory(path);

	FILE *pFile = fopen((dir + "/run.log").c_str(), "a");
	if (pFile == NULL)
		return false;

	get_time(start, currentTime, hours, minutes, seconds);
	fprintf(pFile, ("%4.0fh:%2.0fm:%2.0fs %s\n"), hours, minutes, seconds, text.c_str());
	fclose(pFile);
	if (verbose)
		fdtd::io::log("run", "%s, elapsed %.0f h %.0f m %.0f s", text.c_str(), hours, minutes, seconds);

	return error;
}

// STRUCTURED FIELD OUTPUT (HDF5)

namespace
{
const char* const kAxis = "XYZ";

std::vector<int> domain_size(int const* globalInnerSize)
{
	return {globalInnerSize[0], globalInnerSize[1], globalInnerSize[2]};
}
} // namespace

/**
 * @brief Writes the three components of a 3D field on one plane to an HDF5 file.
 *
 * The file holds datasets /x, /y, /z of shape [n_B][n_A] (A < B, e.g. XZ -> [nz][nx]), x fastest,
 * plus attributes domain_size = [nx, ny, nz], physical domain, no PMLs, plane, and position, the
 * plane's cell along its normal.
 * Written collectively by the ranks of outputCommPlotPlane.
 *
 * @param filename Output path (".h5")
 * @param arr_x, arr_y, arr_z Components of the 3D field
 * @param outputCommPlotPlane Ranks whose chunks intersect the plane and are not pure CPML
 * @param normalToPlaneIdx Index of the plane along its normal
 * @param startindices Global, physical-domain, index of this rank's first cell
 * @param globalInnerSize Cells of the physical domain in x, y, z
 * @param nrInnerThis Cells of this rank inside the physical domain
 * @param nrShift Offset from the CPML-padded to the physical domain
 * @param idx_coord_A, idx_coord_B Axes spanning the plane (0 = X, 1 = Y, 2 = Z), A < B
 * @param use_global_coords Arrays are indexed globally, physical domain
 * @param ijk_order Arrays are arr[i][j][k], false: arr[k][j][i]
 * @param with_PMLs Arrays include the CPML padding
 */
template <typename T>
void save_as_binary_2d(string const &filename, T ***arr_x, T ***arr_y, T ***arr_z, MPI_Comm const &outputCommPlotPlane, int const &normalToPlaneIdx,
int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift, int const &idx_coord_A, int const &idx_coord_B,
bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs)
{
	int idxShift[3] = {0, 0, 0};
	if (use_global_coords)
		for (int d = 0; d < 3; ++d)
			idxShift[d] = startindices[d];
	if (with_PMLs)
		for (int d = 0; d < 3; ++d)
			idxShift[d] += nrShift[d];

	const int na = nrInnerThis[idx_coord_A];
	const int nb = nrInnerThis[idx_coord_B];
	const int normal = 3 - idx_coord_A - idx_coord_B;
	std::vector<T> plane_x(size_t(na) * nb);
	std::vector<T> plane_y(plane_x.size());
	std::vector<T> plane_z(plane_x.size());

	for (int a = 0; a < na; ++a)
		for (int b = 0; b < nb; ++b)
		{
			int g[3];
			g[idx_coord_A] = a + idxShift[idx_coord_A];
			g[idx_coord_B] = b + idxShift[idx_coord_B];
			g[normal] = normalToPlaneIdx + idxShift[normal];
			const size_t n = size_t(a) * nb + b;
			if (ijk_order)
			{
				plane_x[n] = arr_x[g[0]][g[1]][g[2]];
				plane_y[n] = arr_y[g[0]][g[1]][g[2]];
				plane_z[n] = arr_z[g[0]][g[1]][g[2]];
			}
			else
			{
				plane_x[n] = arr_x[g[2]][g[1]][g[0]];
				plane_y[n] = arr_y[g[2]][g[1]][g[0]];
				plane_z[n] = arr_z[g[2]][g[1]][g[0]];
			}
		}

	const auto block = h5io::Block::plane(globalInnerSize[idx_coord_A], globalInnerSize[idx_coord_B],
										  startindices[idx_coord_A], startindices[idx_coord_B], na, nb);
	h5io::File file(filename, outputCommPlotPlane);
	file.attribute("domain_size", domain_size(globalInnerSize));
	file.attribute("plane", std::string{kAxis[idx_coord_A], kAxis[idx_coord_B]});
	file.attribute("position", std::vector<int>{normalToPlaneIdx + startindices[normal]});
	h5io::write(file, "x", block, plane_x.data());
	h5io::write(file, "y", block, plane_y.data());
	h5io::write(file, "z", block, plane_z.data());
}

template void save_as_binary_2d<int>(string const &filename, int ***arr_x, int ***arr_y, int ***arr_z, MPI_Comm const &outputCommPlotPlane,
int const &normalToPlaneIdx, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift,
int const &idx_coord_A, int const &idx_coord_B, bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);

template void save_as_binary_2d<double>(string const &filename, double ***arr_x, double ***arr_y, double ***arr_z, MPI_Comm const &outputCommPlotPlane,
int const &normalToPlaneIdx, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift,
int const &idx_coord_A, int const &idx_coord_B, bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);

template void save_as_binary_2d<std::complex<double>>(string const &filename, std::complex<double> ***arr_x, std::complex<double> ***arr_y, std::complex<double> ***arr_z, MPI_Comm const &outputCommPlotPlane,
int const &normalToPlaneIdx, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift,
int const &idx_coord_A, int const &idx_coord_B, bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);

/**
 * @brief Writes the three components of a 3D field to an HDF5 file.
 *
 * Datasets /x, /y, /z of shape [nz][ny][nx], x fastest, physical domain, no PMLs, plus the
 * attribute domain_size. Written collectively by the ranks of outputCommPlot.
 * Parameters as for save_as_binary_2d.
 */
template <typename T>
void save_as_binary_3d(string const &filename, T ***arr_x, T ***arr_y, T ***arr_z,  MPI_Comm const &outputCommPlot,
int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift,
bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs)
{
	int idxShift[3] = {0, 0, 0};
	if (use_global_coords)
		for (int d = 0; d < 3; ++d)
			idxShift[d] = startindices[d];
	if (with_PMLs)
		for (int d = 0; d < 3; ++d)
			idxShift[d] += nrShift[d];

	const int nx = nrInnerThis[0];
	const int ny = nrInnerThis[1];
	const int nz = nrInnerThis[2];
	std::vector<T> vol_x(size_t(nx) * ny * nz);
	std::vector<T> vol_y(vol_x.size());
	std::vector<T> vol_z(vol_x.size());
	for (int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
			for (int k = 0; k < nz; ++k)
			{
				const size_t n = (size_t(i) * ny + j) * nz + k;
				const int I = i + idxShift[0];
				const int J = j + idxShift[1];
				const int K = k + idxShift[2];
				if (ijk_order)
				{
					vol_x[n] = arr_x[I][J][K];
					vol_y[n] = arr_y[I][J][K];
					vol_z[n] = arr_z[I][J][K];
				}
				else
				{
					vol_x[n] = arr_x[K][J][I];
					vol_y[n] = arr_y[K][J][I];
					vol_z[n] = arr_z[K][J][I];
				}
			}

	const hsize_t global[3] = {hsize_t(globalInnerSize[0]), hsize_t(globalInnerSize[1]), hsize_t(globalInnerSize[2])};
	const hsize_t offset[3] = {hsize_t(startindices[0]), hsize_t(startindices[1]), hsize_t(startindices[2])};
	const hsize_t count[3] = {hsize_t(nx), hsize_t(ny), hsize_t(nz)};
	const auto block = h5io::Block::volume(global, offset, count);
	h5io::File file(filename, outputCommPlot);
	file.attribute("domain_size", domain_size(globalInnerSize));
	h5io::write(file, "x", block, vol_x.data());
	h5io::write(file, "y", block, vol_y.data());
	h5io::write(file, "z", block, vol_z.data());
}

template void save_as_binary_3d<int>(string const &filename, int ***arr_x, int ***arr_y, int ***arr_z,
MPI_Comm const &outputComm, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift,
bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);

template void save_as_binary_3d<double>(string const &filename, double ***arr_x, double ***arr_y, double ***arr_z,
MPI_Comm const &outputComm, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift,
bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);

template void save_as_binary_3d<std::complex<double>>(string const &filename, std::complex<double> ***arr_x, std::complex<double> ***arr_y, std::complex<double> ***arr_z,
MPI_Comm const &outputComm, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift,
bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);

void fdtd::io::log(const char *key, const char *format, ...)
{
	if (log_muted)
		return;
	int rank = 0;
	int mpi_on = 0;
	int mpi_off = 0; // this process, 0 outside MPI
	MPI_Initialized(&mpi_on);
	MPI_Finalized(&mpi_off);
	if (mpi_on && !mpi_off)
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	char stack[1024];
	std::va_list args;
	va_start(args, format);
	std::va_list copy;
	va_copy(copy, args);
	const int n = std::vsnprintf(stack, sizeof stack, format, args);
	va_end(args);
	std::string text;
	if (n >= int(sizeof stack))
	{
		text.resize(n + 1);
		std::vsnprintf(text.data(), n + 1, format, copy);
		text.resize(n);
	}
	else if (n > 0)
		text.assign(stack, n);
	va_end(copy);

	std::string out;
	std::size_t start = 0;
	while (start <= text.size())
	{
		std::size_t end = text.find('\n', start);
		if (end == std::string::npos)
			end = text.size();
		std::string line = text.substr(start, end - start);
		start = end + 1;
		for (char &ch : line)
			if (ch == '\t' || ch == '\r')
				ch = ' ';
		const std::size_t first = line.find_first_not_of(' ');
		const std::size_t last = line.find_last_not_of(' ');
		if (first == std::string::npos)
			continue;
		line = line.substr(first, last - first + 1);
		if (line.find_first_not_of("-=*#_~") == std::string::npos)
			continue; // a rule
		out += key;
		if (rank != 0)
			out += "[" + std::to_string(rank) + "]";
		out += ": " + line + "\n";
	}
	std::fputs(out.c_str(), stdout);
}
