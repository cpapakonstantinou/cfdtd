/**
 * @file io/h5io.cc
 * @brief Parallel HDF5 writers, see h5io.h.
 */
#include <algorithm>
#include <type_traits>
#include "io/h5io.h"

#include <cstdio>
#include <stdexcept>

namespace h5io
{
namespace
{

/// HDF5 failures inside a collective write cannot be recovered per rank:
/// report and abort the whole job instead of deadlocking the others.
[[noreturn]] void fail(const std::string& what, const std::string& path)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::fprintf(stderr, "error: rank %d: HDF5: %s: %s\n", rank, path.c_str(), what.c_str());
	H5Eprint2(H5E_DEFAULT, stderr);
	MPI_Abort(MPI_COMM_WORLD, 3);
	std::abort();
}

void check(herr_t status, const char* what, const std::string& path)
{
	if (status < 0)
		fail(what, path);
}

hid_t check_id(hid_t id, const char* what, const std::string& path)
{
	if (id < 0)
		fail(what, path);
	return id;
}

template <typename T> struct type_of;
template <> struct type_of<int>    { static hid_t mem() { return H5Tcopy(H5T_NATIVE_INT); }    static hid_t file() { return H5Tcopy(H5T_STD_I32LE); } };
template <> struct type_of<float>  { static hid_t mem() { return H5Tcopy(H5T_NATIVE_FLOAT); }  static hid_t file() { return H5Tcopy(H5T_IEEE_F32LE); } };
template <> struct type_of<double> { static hid_t mem() { return H5Tcopy(H5T_NATIVE_DOUBLE); } static hid_t file() { return H5Tcopy(H5T_IEEE_F64LE); } };

} // namespace

namespace
{
/**
 * @brief Visits this rank's block in file order: the solver's first axis varies fastest in the file
 *
 * f(source, destination): the offset of an element in the buffer laid out as block.local (row-major,
 * solver order) and its offset in a dense buffer of the block in file order. Blocks of rank 1 or 2
 * are padded in front to rank 3. Tiles keep the strided side in cache.
 */
template <typename F>
void for_each_reordered(const Block &b, F f)
{
	const size_t rank = b.global.size();
	const std::vector<hsize_t> &local = b.local.empty() ? b.count : b.local;
	hsize_t c[3] = {1, 1, 1};
	hsize_t l[3] = {1, 1, 1};
	hsize_t s[3] = {0, 0, 0};
	for (size_t d = 0; d < rank; ++d)
	{
		c[3 - rank + d] = b.count[d];
		l[3 - rank + d] = local[d];
		s[3 - rank + d] = b.start.empty() ? 0 : b.start[d];
	}
	const hsize_t tile = 32;
	for (hsize_t i1 = 0; i1 < c[1]; ++i1)
		for (hsize_t t0 = 0; t0 < c[0]; t0 += tile)
			for (hsize_t t2 = 0; t2 < c[2]; t2 += tile)
				for (hsize_t i0 = t0; i0 < std::min(t0 + tile, c[0]); ++i0)
				{
					const hsize_t row = ((s[0] + i0) * l[1] + s[1] + i1) * l[2] + s[2];
					for (hsize_t i2 = t2; i2 < std::min(t2 + tile, c[2]); ++i2)
						f(row + i2, (i2 * c[1] + i1) * c[0] + i0);
				}
}

std::vector<hsize_t> reversed(const std::vector<hsize_t> &v)
{
	return std::vector<hsize_t>(v.rbegin(), v.rend());
}
} // namespace

/// Shared dataset path, friend of File for the transfer property list. Buffers are dense, in file order.
struct Writer
{
	template <typename T>
	static void write(File &f, const std::string &name, const Block &b, const T *dense)
	{
		const int rank = int(b.global.size());
		const std::string where = f.path() + ":" + name;
		const std::vector<hsize_t> global = reversed(b.global);
		const std::vector<hsize_t> offset = reversed(b.offset);
		const std::vector<hsize_t> count = reversed(b.count);

		hid_t ftype = type_of<T>::file();
		hid_t mtype = type_of<T>::mem();
		hid_t fspace = check_id(H5Screate_simple(rank, global.data(), nullptr), "H5Screate_simple", where);
		hid_t dset = check_id(H5Dcreate2(f._file, name.c_str(), ftype, fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), "H5Dcreate2", where);
		hid_t mspace = check_id(H5Screate_simple(rank, count.data(), nullptr), "H5Screate_simple", where);
		if (b.elements() == 0)
		{
			// ranks without data still take part in the collective write
			H5Sselect_none(fspace);
			H5Sselect_none(mspace);
		}
		else
			check(H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset.data(), nullptr, count.data(), nullptr), "file hyperslab", where);
		check(H5Dwrite(dset, mtype, mspace, fspace, f._xfer, dense), "H5Dwrite", where);

		H5Sclose(mspace);
		H5Dclose(dset);
		H5Sclose(fspace);
		H5Tclose(mtype);
		H5Tclose(ftype);
		f._datasets.push_back({name, global, std::is_same<T, int>::value ? "Int" : "Float", int(sizeof(T))});
	}

	template <typename T>
	static void read(File &f, const std::string &name, const Block &b, T *dense)
	{
		const int rank = int(b.global.size());
		const std::string where = f.path() + ":" + name;
		const std::vector<hsize_t> global = reversed(b.global);
		const std::vector<hsize_t> offset = reversed(b.offset);
		const std::vector<hsize_t> count = reversed(b.count);
		hid_t dset = check_id(H5Dopen2(f._file, name.c_str(), H5P_DEFAULT), "H5Dopen2", where);
		hid_t fspace = H5Dget_space(dset);
		std::vector<hsize_t> extent(rank, 0);
		if (H5Sget_simple_extent_ndims(fspace) != rank || H5Sget_simple_extent_dims(fspace, extent.data(), nullptr) < 0 || extent != global)
			fail("dataset extent does not match the simulation domain", where);

		hid_t mtype = type_of<T>::mem();
		hid_t mspace = check_id(H5Screate_simple(rank, count.data(), nullptr), "H5Screate_simple", where);
		if (b.elements() == 0)
		{
			H5Sselect_none(fspace);
			H5Sselect_none(mspace);
		}
		else
			check(H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset.data(), nullptr, count.data(), nullptr), "file hyperslab", where);
		check(H5Dread(dset, mtype, mspace, fspace, f._xfer, dense), "H5Dread", where);
		H5Sclose(mspace);
		H5Tclose(mtype);
		H5Sclose(fspace);
		H5Dclose(dset);
	}
};

// Block

Block Block::plane(hsize_t ga, hsize_t gb, hsize_t oa, hsize_t ob, hsize_t ca, hsize_t cb)
{
	return Block{{ga, gb}, {oa, ob}, {ca, cb}, {}, {}};
}

Block Block::plane(hsize_t ga, hsize_t gb, hsize_t oa, hsize_t ob, hsize_t ca, hsize_t cb,
				   hsize_t rows, hsize_t row_length, hsize_t sa, hsize_t sb)
{
	return Block{{ga, gb}, {oa, ob}, {ca, cb}, {rows, row_length}, {sa, sb}};
}

Block Block::volume(const hsize_t global[3], const hsize_t offset[3], const hsize_t count[3])
{
	return Block{{global[0], global[1], global[2]},
				 {offset[0], offset[1], offset[2]},
				 {count[0], count[1], count[2]},
				 {},
				 {}};
}

hsize_t Block::elements() const
{
	hsize_t n = 1;
	for (hsize_t c : count)
		n *= c;
	return n;
}

// File

namespace
{
struct XdmfSettings
{
	bool on = false;
	double cell[3] = {1, 1, 1};
} xdmf;
} // namespace

void set_xdmf(bool on, double dx, double dy, double dz)
{
	xdmf.on = on;
	xdmf.cell[0] = dx;
	xdmf.cell[1] = dy;
	xdmf.cell[2] = dz;
}

File::File(const std::string& path, MPI_Comm comm, Access access)
	: _path(path), _comm(comm), _writing(access == Access::create)
{
	H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr); // errors are reported by fail()
	hid_t fapl = check_id(H5Pcreate(H5P_FILE_ACCESS), "H5Pcreate", path);
	check(H5Pset_fapl_mpio(fapl, comm, MPI_INFO_NULL), "H5Pset_fapl_mpio", path);
	_file = access == Access::create
				? check_id(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl), "H5Fcreate", path)
				: check_id(H5Fopen(path.c_str(), H5F_ACC_RDONLY, fapl), "H5Fopen", path);
	H5Pclose(fapl);
	_xfer = check_id(H5Pcreate(H5P_DATASET_XFER), "H5Pcreate", path);
	check(H5Pset_dxpl_mpio(_xfer, H5FD_MPIO_COLLECTIVE), "H5Pset_dxpl_mpio", path);
}

File::~File()
{
	if (_xfer >= 0)
		H5Pclose(_xfer);
	if (_file >= 0)
		H5Fclose(_file);
	int rank = 0;
	if (_writing && xdmf.on && _comm != MPI_COMM_NULL && MPI_Comm_rank(_comm, &rank) == MPI_SUCCESS && rank == 0)
		write_xdmf();
}

void File::write_xdmf() const
{
	// the solver axis of each dimension of a dataset as stored, slowest first: the reverse of the solver order
	auto axes = [this](size_t rank, int out[3]) {
		if (rank == 3)
		{
			out[0] = 2, out[1] = 1, out[2] = 0;
			return true;
		}
		if (rank != 2)
			return false;
		const std::string p = _plane;
		if (p.size() == 2)
		{
			out[0] = p[1] - 'X', out[1] = p[0] - 'X';
			return out[0] >= 0 && out[0] < 3 && out[1] >= 0 && out[1] < 3;
		}
		const int n = _face.empty() ? 2 : _face[0] - 'x';
		out[0] = n == 2 ? 1 : 2, out[1] = n == 0 ? 1 : 0;
		return true;
	};
	std::vector<std::vector<hsize_t>> grids; // distinct extents, in order of first use
	for (const Dataset& d : _datasets)
	{
		int a[3] = {0, 1, 2};
		if (axes(d.extent.size(), a) && std::find(grids.begin(), grids.end(), d.extent) == grids.end())
			grids.push_back(d.extent);
	}
	if (grids.empty())
		return;

	const std::string xmf = (_path.size() > 3 && _path.compare(_path.size() - 3, 3, ".h5") == 0 ? _path.substr(0, _path.size() - 3) : _path) + ".xmf";
	const std::string h5 = _path.substr(_path.find_last_of('/') + 1);
	FILE* f = std::fopen(xmf.c_str(), "w");
	if (!f)
		return;
	std::fprintf(f, "<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
	std::fprintf(f, "<!-- %s: viewer X, Y, Z = solver x, y, z; arrays stored [z][y][x], x fastest -->\n", h5.c_str());
	std::fprintf(f, "<Xdmf Version=\"2.0\">\n <Domain>\n");
	for (size_t g = 0; g < grids.size(); ++g)
	{
		const std::vector<hsize_t>& e = grids[g];
		int a[3] = {0, 1, 2};
		axes(e.size(), a);
		// a 3D grid, z y x slowest first; a plane is one point thick along its normal, at its position
		hsize_t points[3] = {1, 1, 1};
		double origin[3] = {0, 0, 0};
		for (size_t d = 0; d < e.size(); ++d)
			points[a[d]] = e[d];
		if (e.size() == 2)
		{
			const int normal = 3 - a[0] - a[1];
			origin[normal] = _position * xdmf.cell[normal];
		}
		const std::string dims = std::to_string(points[2]) + " " + std::to_string(points[1]) + " " + std::to_string(points[0]);
		std::fprintf(f, "  <Grid Name=\"grid_%zu\" GridType=\"Uniform\">\n", g);
		std::fprintf(f, "   <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%s\"/>\n", dims.c_str());
		std::fprintf(f, "   <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
		std::fprintf(f, "    <DataItem Dimensions=\"3\" Format=\"XML\">%.9g %.9g %.9g</DataItem>\n", origin[2], origin[1], origin[0]);
		std::fprintf(f, "    <DataItem Dimensions=\"3\" Format=\"XML\">%.9g %.9g %.9g</DataItem>\n", xdmf.cell[2], xdmf.cell[1], xdmf.cell[0]);
		std::fprintf(f, "   </Geometry>\n");
		for (const Dataset& d : _datasets)
			if (d.extent == e)
			{
				std::fprintf(f, "   <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", d.name.c_str());
				std::fprintf(f, "    <DataItem Dimensions=\"%s\" NumberType=\"%s\" Precision=\"%d\" Format=\"HDF\">%s:/%s</DataItem>\n",
							 dims.c_str(), d.number, d.precision, h5.c_str(), d.name.c_str());
				std::fprintf(f, "   </Attribute>\n");
			}
		std::fprintf(f, "  </Grid>\n");
	}
	std::fprintf(f, " </Domain>\n</Xdmf>\n");
	std::fclose(f);
}

void File::attribute(const std::string& name, const std::vector<int>& values)
{
	if (name == "position" && values.size() == 1)
		_position = values[0];
	hsize_t n = values.size();
	hid_t space = H5Screate_simple(1, &n, nullptr);
	hid_t attr = check_id(H5Acreate2(_file, name.c_str(), H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT),
						  "H5Acreate2", _path + ":" + name);
	check(H5Awrite(attr, H5T_NATIVE_INT, values.data()), "H5Awrite", _path + ":" + name);
	H5Aclose(attr);
	H5Sclose(space);
}

void File::attribute(const std::string& name, const std::string& value)
{
	if (name == "plane")
		_plane = value;
	if (name == "face")
		_face = value;
	hid_t type = H5Tcopy(H5T_C_S1);
	H5Tset_size(type, value.size() ? value.size() : 1);
	hid_t space = H5Screate(H5S_SCALAR);
	hid_t attr = check_id(H5Acreate2(_file, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT),
						  "H5Acreate2", _path + ":" + name);
	check(H5Awrite(attr, type, value.c_str()), "H5Awrite", _path + ":" + name);
	H5Aclose(attr);
	H5Sclose(space);
	H5Tclose(type);
}

void File::attribute(const std::string& name, double value)
{
	hid_t space = H5Screate(H5S_SCALAR);
	hid_t attr = check_id(H5Acreate2(_file, name.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT),
						  "H5Acreate2", _path + ":" + name);
	check(H5Awrite(attr, H5T_NATIVE_DOUBLE, &value), "H5Awrite", _path + ":" + name);
	H5Aclose(attr);
	H5Sclose(space);
}

std::vector<int> File::int_attribute(const std::string& name)
{
	const std::string where = _path + ":" + name;
	hid_t attr = check_id(H5Aopen(_file, name.c_str(), H5P_DEFAULT), "H5Aopen", where);
	hid_t space = H5Aget_space(attr);
	std::vector<int> v(size_t(H5Sget_simple_extent_npoints(space)));
	check(H5Aread(attr, H5T_NATIVE_INT, v.data()), "H5Aread", where);
	H5Sclose(space);
	H5Aclose(attr);
	return v;
}

// writers

template <typename T>
void write(File &file, const std::string &name, const Block &block, const T *data)
{
	std::vector<T> dense(block.elements());
	for_each_reordered(block, [&](hsize_t from, hsize_t to) { dense[to] = data[from]; });
	Writer::write(file, name, block, dense.data());
}

template <typename T>
void write_complex(File &file, const std::string &name, const Block &block, const T *re, const T *im)
{
	write(file, name + "_re", block, re);
	write(file, name + "_im", block, im);
}

/// std::complex buffers: real and imaginary parts as two datasets
template <>
void write<std::complex<double>>(File &file, const std::string &name, const Block &block, const std::complex<double> *data)
{
	std::vector<double> re(block.elements());
	std::vector<double> im(block.elements());
	for_each_reordered(block, [&](hsize_t from, hsize_t to) {
		re[to] = data[from].real();
		im[to] = data[from].imag();
	});
	Writer::write(file, name + "_re", block, re.data());
	Writer::write(file, name + "_im", block, im.data());
}

template <typename T>
void read(File &file, const std::string &name, const Block &block, T *data)
{
	std::vector<T> dense(block.elements());
	Writer::read(file, name, block, dense.data());
	for_each_reordered(block, [&](hsize_t to, hsize_t from) { data[to] = dense[from]; });
}

template void read<double>(File&, const std::string&, const Block&, double*);
template void read<float>(File&, const std::string&, const Block&, float*);
template void write<int>(File&, const std::string&, const Block&, const int*);
template void write<float>(File&, const std::string&, const Block&, const float*);
template void write<double>(File&, const std::string&, const Block&, const double*);
template void write_complex<double>(File&, const std::string&, const Block&, const double*, const double*);
template void write_complex<float>(File&, const std::string&, const Block&, const float*, const float*);

} // namespace h5io
