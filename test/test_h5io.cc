// Parallel HDF5 writer/reader tests. Run with any number of ranks:
//   mpirun -np 4 ./test_h5io
#include "io/h5io.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>
#include <vector>

namespace
{
int _rank = 0;
int _size = 1;
int failures = 0;

void expect(bool ok, const char *what)
{
	int all = ok ? 1 : 0;
	int global = 0;
	MPI_Allreduce(&all, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	if (!global)
	{
		if (_rank == 0)
			std::printf("FAIL: %s\n", what);
		++failures;
	}
}

std::string temp_path(const char *name)
{
	// same path on every rank: rank 0's pid keeps concurrent test runs apart
	int pid = int(getpid());
	MPI_Bcast(&pid, 1, MPI_INT, 0, MPI_COMM_WORLD);
	return std::string("/tmp/test_h5io_") + std::to_string(pid) + "_" + name + ".h5";
}

// Global n0 x n1 plane split in rows across ranks; each rank keeps a padded buffer
// with one ghost row/column and a one-cell margin before the data.
void padded_plane_round_trip()
{
	const hsize_t n0 = 7;
	const hsize_t n1 = 5;
	const hsize_t first = n0 * _rank / _size;
	const hsize_t last = n0 * (_rank + 1) / _size;
	const hsize_t rows = last - first; // may be 0 for some ranks when size > n0
	const hsize_t lrows = rows + 2;
	const hsize_t lcols = n1 + 2;
	std::vector<double> local(lrows * lcols, -1.0);
	for (hsize_t a = 0; a < rows; ++a)
		for (hsize_t b = 0; b < n1; ++b)
			local[(a + 1) * lcols + b + 1] = double((first + a) * 100 + b);

	const std::string path = temp_path("plane");
	{
		h5io::File f(path, MPI_COMM_WORLD);
		f.attribute("domain_size", std::vector<int>{int(n0), 1, int(n1)});
		f.attribute("plane", std::string("XZ"));
		h5io::write(f, "x", h5io::Block::plane(n0, n1, first, 0, rows, n1, lrows, lcols, 1, 1), local.data());
	}
	// read back the whole plane on every rank
	std::vector<double> all(n0 * n1, 0.0);
	{
		h5io::File f(path, MPI_COMM_WORLD, h5io::File::Access::read);
		expect(f.int_attribute("domain_size") == std::vector<int>{int(n0), 1, int(n1)}, "plane: domain_size attribute");
		h5io::read(f, "x", h5io::Block::plane(n0, n1, 0, 0, n0, n1), all.data());
	}
	bool ok = true;
	for (hsize_t a = 0; a < n0; ++a)
		for (hsize_t b = 0; b < n1; ++b)
			ok &= all[a * n1 + b] == double(a * 100 + b);
	expect(ok, "plane: padded window round trip");
	if (_rank == 0)
		std::remove(path.c_str());
}

// Complex dataset from separate Re/Im buffers, 3D, split along the first axis.
void complex_volume()
{
	const hsize_t g[3] = {4, 3, 2};
	const hsize_t first = g[0] * _rank / _size;
	const hsize_t last = g[0] * (_rank + 1) / _size;
	const hsize_t off[3] = {first, 0, 0};
	const hsize_t cnt[3] = {last - first, g[1], g[2]};
	std::vector<double> re(cnt[0] * g[1] * g[2]);
	std::vector<double> im(re.size());
	for (hsize_t i = 0; i < cnt[0]; ++i)
		for (hsize_t j = 0; j < g[1]; ++j)
			for (hsize_t k = 0; k < g[2]; ++k)
			{
				const hsize_t n = (i * g[1] + j) * g[2] + k;
				re[n] = double(first + i) + 0.25 * double(j);
				im[n] = -double(k) - 0.5;
			}
	const std::string path = temp_path("complex");
	h5io::set_xdmf(true, 1e-9, 2e-9, 3e-9);
	{
		h5io::File f(path, MPI_COMM_WORLD);
		h5io::write_complex(f, "z", h5io::Block::volume(g, off, cnt), re.data(), im.data());
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (_rank == 0)
	{
		// independent serial read with the plain HDF5 API: z_re and z_im
		hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		std::vector<double> r(g[0] * g[1] * g[2]);
		std::vector<double> q(g[0] * g[1] * g[2]);
		hid_t dr = H5Dopen2(file, "z_re", H5P_DEFAULT);
		hid_t di = H5Dopen2(file, "z_im", H5P_DEFAULT);
		const bool opened = dr >= 0 && di >= 0;
		if (opened)
		{
			H5Dread(dr, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, r.data());
			H5Dread(di, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, q.data());
			H5Dclose(dr);
			H5Dclose(di);
		}
		H5Fclose(file);
		bool ok = opened;
		for (hsize_t i = 0; i < g[0]; ++i)
			for (hsize_t j = 0; j < g[1]; ++j)
				for (hsize_t k = 0; k < g[2]; ++k)
				{
					const hsize_t n = (i * g[1] + j) * g[2] + k;
					ok &= r[n] == double(i) + 0.25 * double(j) && q[n] == -double(k) - 0.5;
				}
		if (!ok)
			std::printf("FAIL: complex volume values\n"), ++failures;
		// the XDMF companion, on by default in this test: set_xdmf below
		const std::string xmf = path.substr(0, path.size() - 3) + ".xmf";
		std::FILE *x = std::fopen(xmf.c_str(), "r");
		std::string text;
		if (x)
		{
			char buf[4096];
			size_t got;
			while ((got = std::fread(buf, 1, sizeof buf, x)) > 0)
				text.append(buf, got);
			std::fclose(x);
		}
		const std::string dims = std::to_string(g[0]) + " " + std::to_string(g[1]) + " " + std::to_string(g[2]);
		if (text.find("3DCoRectMesh") == std::string::npos || text.find("Dimensions=\"" + dims + "\"") == std::string::npos ||
			text.find(":/z_re<") == std::string::npos || text.find(":/z_im<") == std::string::npos)
			std::printf("FAIL: xdmf companion\n"), ++failures;
		std::remove(xmf.c_str());
		std::remove(path.c_str());
	}
	MPI_Bcast(&failures, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

// Integer maps, material indices, keep their type.
void integer_plane()
{
	const int v = _rank + 1;
	const std::string path = temp_path("int");
	{
		h5io::File f(path, MPI_COMM_WORLD);
		h5io::write(f, "x", h5io::Block::plane(_size, 1, _rank, 0, 1, 1), &v);
	}
	std::vector<int> back(_size, 0);
	{
		h5io::File f(path, MPI_COMM_WORLD, h5io::File::Access::read);
		hid_t d = H5Dopen2(f.id(), "x", H5P_DEFAULT);
		hid_t t = H5Dget_type(d);
		expect(H5Tget_class(t) == H5T_INTEGER, "int plane: integer dataset type");
		H5Tclose(t);
		H5Dclose(d);
	}
	if (_rank == 0)
		std::remove(path.c_str());
}
} // namespace

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &_size);
	padded_plane_round_trip();
	complex_volume();
	integer_plane();
	if (_rank == 0)
		std::printf(failures ? "%d h5io test(s) failed\n" : "all h5io tests passed (%d ranks)\n", failures ? failures : _size);
	MPI_Finalize();
	return failures ? 1 : 0;
}
