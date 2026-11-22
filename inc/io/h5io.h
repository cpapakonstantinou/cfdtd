/**
 * @file io/h5io.h
 * @brief Parallel HDF5 output for structured, plane / volume, data.
 *
 * Every structured output is one HDF5 file written collectively by the ranks
 * of an MPI communicator through the MPI-IO driver. Each rank contributes one
 * rectangular block of every dataset:
 *
 *     global  extent of the dataset in the file           (e.g. nx, ny)
 *     offset  where this rank's block starts in the file
 *     count   size of this rank's block
 *     local   allocated extent of the rank's memory buffer, may include
 *             ghost cells / padding
 *     start   where the block starts inside the memory buffer
 *
 * Blocks are given in the solver's order, [x][y][z] or the plane's [a][b], a < b. The file stores
 * the reverse, x fastest: [z][y][x], [b][a]; write and read reorder through a dense buffer, in
 * tiles, so readers and viewers see x across and the last axis up.
 *
 * Complex data is stored as two real datasets <name>_re and <name>_im, ParaView's
 * XDMF reader cannot read HDF5 compound types.
 *
 * XDMF: with set_xdmf(true, ...) every file written here gets a <file>.xmf next to
 * it, written by rank 0 of the file's communicator when the file closes, which
 * describes the 2D and 3D datasets as attributes on a uniform grid, so the outputs
 * open in ParaView / VisIt. The grid has the cell size given to set_xdmf; datasets
 * of the same extent share one grid. Viewer X, Y, Z = solver x, y, z: every grid is
 * 3D, a plane one point thick along its normal at its "position" attribute (0 when
 * it has none). Verified with VTK's vtkXdmfReader, the ParaView reader: every point
 * at its solver coordinates. 2D datasets use the axes
 * named by the file's "plane" attribute, "XZ": first index x, second z, or the two
 * axes tangential to its "face" attribute; 1D datasets are not described.
 *
 * Layout written by the field writers (save_as_binary_2d/3D):
 *     /x, /y, /z            one dataset per vector component, x_re, x_im, ... if complex
 *     attribute domain_size [nx, ny, nz] of the physical domain, no PMLs
 *     attribute plane       "XY" | "XZ" | "YZ", plane outputs only
 *     attribute position    the plane's cell along its normal, plane outputs only
 */
#ifndef FDTD_H5IO_H
#define FDTD_H5IO_H

#include <complex>
#include <string>
#include <vector>

#include <hdf5.h>
#ifndef H5_HAVE_PARALLEL
#error "this HDF5 has no MPI support: load a parallel HDF5, or name its pkg-config package with make HDF5_PKG=..."
#endif
#include <mpi.h>

namespace h5io
{

/// XDMF companions on or off, and the cell size [m] of the grid they describe.
void set_xdmf(bool on, double dx, double dy, double dz);

/// Placement of one rank's block inside a dataset, see file comment.
struct Block
{
	std::vector<hsize_t> global;
	std::vector<hsize_t> offset;
	std::vector<hsize_t> count;
	std::vector<hsize_t> local; ///< empty: same as count, dense buffer
	std::vector<hsize_t> start; ///< empty: all zeros

	/// Dense 2D block, buffer is exactly count[0] x count[1].
	static Block plane(hsize_t ga, hsize_t gb, hsize_t oa, hsize_t ob, hsize_t ca, hsize_t cb);
	/// 2D window of a padded buffer: rows x row_length elements, block at (sa, sb).
	static Block plane(hsize_t ga, hsize_t gb, hsize_t oa, hsize_t ob, hsize_t ca, hsize_t cb,
					   hsize_t rows, hsize_t row_length, hsize_t sa = 0, hsize_t sb = 0);
	/// Dense 3D block.
	static Block volume(const hsize_t global[3], const hsize_t offset[3], const hsize_t count[3]);

	hsize_t elements() const;
};

/// An HDF5 file opened collectively on `comm`. Closed on destruction.
class File
{
public:
	enum class Access
	{
		create, ///< create, truncating an existing file
		read    ///< open an existing file read-only
	};
	File(const std::string& path, MPI_Comm comm, Access access = Access::create);
	~File();
	File(const File&) = delete;
	File& operator=(const File&) = delete;

	hid_t id() const { return _file; }
	const std::string& path() const { return _path; }

	/// File-level attributes, collective: every rank passes the same value.
	void attribute(const std::string& name, const std::vector<int>& values);
	void attribute(const std::string& name, const std::string& value);
	void attribute(const std::string& name, double value);
	/// Reads an integer vector attribute, collective.
	std::vector<int> int_attribute(const std::string& name);

private:
	struct Dataset
	{
		std::string name;
		std::vector<hsize_t> extent;
		const char* number; ///< XDMF NumberType
		int precision;      ///< bytes
	};
	void write_xdmf() const;

	hid_t _file = -1;
	hid_t _xfer = -1;
	std::string _path;
	MPI_Comm _comm = MPI_COMM_NULL;
	bool _writing = false;
	std::vector<Dataset> _datasets;
	std::string _plane; ///< the "plane" / "face" attributes, 2D axes
	std::string _face; ///< the "plane" / "face" attributes, 2D axes
	int _position = 0; ///< the "position" attribute: cell of a plane along its normal
	std::string _description; ///< the "description" attribute, in the XDMF comment
	friend struct Writer;
};

/// Collectively write `data`, laid out as block.local, row-major, to dataset `name`.
template <typename T>
void write(File& file, const std::string& name, const Block& block, const T* data);

/// Collectively read this rank's block of dataset `name` into `data`, laid out as
/// block.local. The dataset extent must equal block.global.
template <typename T>
void read(File& file, const std::string& name, const Block& block, T* data);

/// Collectively write the complex dataset `name`, as <name>_re and <name>_im, from
/// separate real and imaginary buffers, both laid out as block.local.
template <typename T>
void write_complex(File& file, const std::string& name, const Block& block, const T* re, const T* im);

/**
 * What a frequency-domain field file holds: the stored DFT times `scale`, and the attributes that say
 * so. `description` also goes into the comment of the XDMF file.
 *
 *     attribute wavelength                  of the field [m]
 *     attribute normalization               "incident_field" or "pump_field"
 *     attribute harmonic_order              k, the field's frequency over its pump point's
 *     attribute reference_wavelength        of the pump point whose |E_inc| divides [m]
 *     attribute reference_amplitude_v_s_m   that |E_inc| [V s/m]
 *     attribute reference_phase_rad         the phase taken off
 *     attribute dt_s                        the time step: the raw DFT sum is the spectral amplitude over dt
 *     attribute description                 one line
 *
 * The spectral amplitude in V s/m is value * reference_amplitude_v_s_m * exp(i reference_phase_rad).
 */
struct FieldNote
{
	std::complex<double> scale = 1;
	double wavelength = 0;
	std::string normalization;
	int order = 1;
	double reference_wavelength = 0;
	double reference_amplitude = 0;
	double reference_phase = 0;
	double dt = 0;
	std::string description;

	/// collective, the attributes above
	void write(File& file) const;
};

} // namespace h5io

#endif // FDTD_H5IO_H
