/**
 * @file sources/mode_source.cc
 * @brief A waveguide mode from MPB, injected one way on a plane, sources/mode_source.h
 */
#include "sources/mode_source.h"

#include "core/constants.h"
#include "grid/subdomain.h"
#include "io/text.h"
#include "materials/coefficients.h"
#include "materials/models.h"
#include "sources/signal.h"
#include "update/state.h"

#include <hdf5.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>
#include "sources/aux_grid.h"
#include "grid/state.h"
#include "grid/lattice.h"
#include "sources/spec.h"
using fdtd::constants::ROOTPROCESS;
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::fields::Hx;
using fdtd::fields::Hy;
using fdtd::fields::Hz;
using fdtd::grid::subdomain;
using fdtd::materials::DispMateOffsetD;
using fdtd::materials::Drude;
using fdtd::materials::DrudeCP;
using fdtd::materials::Lorentz;
using fdtd::materials::NRMAXMATERIALS;
using fdtd::materials::matind_x;
using fdtd::materials::matind_x_full;
using fdtd::materials::matind_y;
using fdtd::materials::matind_y_full;
using fdtd::materials::matind_z;
using fdtd::materials::matind_z_full;
using fdtd::materials::background_y;

fdtd::sources::ModeSpec fdtd::sources::mode_input;

namespace
{
const fdtd::sources::Spec &sp = fdtd::sources::spec();          ///< what the sources need from the run
const fdtd::grid::Lattice &lat = fdtd::grid::lattice();          ///< the global lattice

/// One complex 2D component of an MPB field file, data[ix * ny + iy]
struct Component
{
	std::vector<std::complex<double>> data;
};

/// The transverse field of one MPB file on MPB's grid: x_i = -Lx/2 + i Lx/nx, likewise y
struct MpbField
{
	int nx = 0;
	int ny = 0;
	double lx = 0; ///< cell size [m]
	double ly = 0; ///< cell size [m]
	Component x;
	Component y;

	/// bilinear interpolation of component c at (px, py) [m] from the cell centre; 0 outside the cell
	double at(const std::vector<double> &c, double px, double py) const
	{
		const double fx = (px + 0.5 * lx) / lx * nx;
		const double fy = (py + 0.5 * ly) / ly * ny;
		const int ix = int(std::floor(fx));
		const int iy = int(std::floor(fy));
		if (ix < 0 || iy < 0 || ix >= nx || iy >= ny)
			return 0;
		const int jx = std::min(ix + 1, nx - 1);
		const int jy = std::min(iy + 1, ny - 1);
		const double ax = fx - ix;
		const double ay = fy - iy;
		return (1 - ax) * ((1 - ay) * c[ix * ny + iy] + ay * c[ix * ny + jy]) + ax * ((1 - ay) * c[jx * ny + iy] + ay * c[jx * ny + jy]);
	}
};

std::vector<double> read_dataset(hid_t file, const std::string &name, std::vector<hsize_t> &dims, const std::string &path)
{
	if (H5Lexists(file, name.c_str(), H5P_DEFAULT) <= 0)
		throw std::runtime_error("mode source: " + path + " has no dataset '" + name + "' (expected MPB field output: x.r, x.i, y.r, ...)");
	const hid_t set = H5Dopen2(file, name.c_str(), H5P_DEFAULT);
	const hid_t space = H5Dget_space(set);
	dims.assign(H5Sget_simple_extent_ndims(space), 0);
	H5Sget_simple_extent_dims(space, dims.data(), nullptr);
	hsize_t n = 1;
	for (hsize_t d : dims)
		n *= d;
	std::vector<double> v(n);
	H5Dread(set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
	H5Sclose(space);
	H5Dclose(set);
	return v;
}

/// Reads the x and y components of an MPB field file, 2D cross-section, and its cell size
MpbField read_mpb(const std::string &path, double unit)
{
	H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
	const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file < 0)
		throw std::runtime_error("mode source: cannot open " + path);
	MpbField f;
	std::vector<hsize_t> dims;
	auto component = [&](const char *axis, Component &c) {
		std::vector<hsize_t> di;
		const auto re = read_dataset(file, std::string(axis) + ".r", dims, path);
		const auto im = read_dataset(file, std::string(axis) + ".i", di, path);
		if (di != dims || re.size() != im.size())
			throw std::runtime_error("mode source: " + path + ": real and imaginary parts differ in size");
		c.data.resize(re.size());
		for (size_t n = 0; n < re.size(); ++n)
			c.data[n] = {re[n], im[n]};
	};
	component("x", f.x);
	component("y", f.y);
	if (dims.size() < 2 || (dims.size() == 3 && dims[2] != 1) || dims.size() > 3)
		throw std::runtime_error("mode source: " + path + " is not a 2D cross-section (solve the mode in MPB with a cell of no size along the waveguide)");
	f.nx = int(dims[0]), f.ny = int(dims[1]);
	// MPB writes the lattice vectors, in units of a, with every field file
	std::vector<hsize_t> ld;
	const auto lattice = read_dataset(file, "lattice vectors", ld, path);
	if (lattice.size() != 9)
		throw std::runtime_error("mode source: " + path + ": 'lattice vectors' is not 3 x 3");
	f.lx = std::hypot(lattice[0], lattice[1], lattice[2]) * unit;
	f.ly = std::hypot(lattice[3], lattice[4], lattice[5]) * unit;
	H5Fclose(file);
	return f;
}

/// Turns the field real: removes the phase of its largest transverse sample
void make_real(MpbField &f, std::vector<double> &x, std::vector<double> &y)
{
	size_t best = 0;
	double top = -1;
	for (size_t n = 0; n < f.x.data.size(); ++n)
	{
		const double a = std::max(std::norm(f.x.data[n]), std::norm(f.y.data[n]));
		if (a > top)
			top = a, best = n;
	}
	const std::complex<double> ref = std::norm(f.x.data[best]) >= std::norm(f.y.data[best]) ? f.x.data[best] : f.y.data[best];
	const std::complex<double> turn = std::abs(ref) > 0 ? std::conj(ref) / std::abs(ref) : 1.0;
	x.resize(f.x.data.size());
	y.resize(f.y.data.size());
	for (size_t n = 0; n < x.size(); ++n)
	{
		x[n] = (f.x.data[n] * turn).real();
		y[n] = (f.y.data[n] * turn).real();
	}
}

/// E update coefficient dt / eps of a material code; false where no plain coefficient applies
bool e_factor(int m, double &c)
{
	const auto &mp = fdtd::materials::model_parameters();
	if (m >= DispMateOffsetD && m < NRMAXMATERIALS)
	{
		const int mm = m - DispMateOffsetD;
		if (mm >= DrudeCP && mm < Lorentz)
		{
			c = mm < Drude ? fdtd::materials::coeffs<fdtd::materials::models::Drude>(mm).ce2 : fdtd::materials::coeffs<fdtd::materials::models::Lorentz>(mm).ce2;
			return true;
		}
		return false;
	}
	c = mp.ce_2[m];
	return true;
}

struct Plane
{
	int n = 1; ///< propagation axis and the transverse axes, cyclic
	int u = 2; ///< propagation axis and the transverse axes, cyclic
	int v = 0; ///< propagation axis and the transverse axes, cyclic
	int e_local = -1; ///< local index of the E plane and of the H plane (p - 1/2), -1: not on this rank
	int h_local = -1; ///< local index of the E plane and of the H plane (p - 1/2), -1: not on this rank
	int nu = 0; ///< local extent of the plane, cells + 1
	int nv = 0; ///< local extent of the plane, cells + 1
	std::vector<double> eu; ///< E_u, E_v of the mode [V/m per unit signal]
	std::vector<double> ev; ///< E_u, E_v of the mode [V/m per unit signal]
	std::vector<double> hu_c; ///< dt / (eps dn) * H_u, H_v: the E corrections
	std::vector<double> hv_c; ///< dt / (eps dn) * H_u, H_v: the E corrections
	std::vector<double> eu_c; ///< dt / (mu0 dn) * E_u, E_v: the H corrections
	std::vector<double> ev_c; ///< dt / (mu0 dn) * E_u, E_v: the H corrections
	double delay = 0;                ///< [s] half a cell at n_eff: the H plane leads the E plane
	bool active = false;
} plane;

std::array<ctensor<double> *, 3> e_fields() { return {&Ex, &Ey, &Ez}; }
std::array<ctensor<double> *, 3> h_fields() { return {&Hx, &Hy, &Hz}; }

/// element (a_n, a_u, a_v) of a field, in the axis order of the plane
double &cell(ctensor<double> &f, int an, int au, int av)
{
	std::array<int, 3> ijk{};
	ijk[plane.n] = an, ijk[plane.u] = au, ijk[plane.v] = av;
	return f[ijk[0]][ijk[1]][ijk[2]];
}

int axis_index(const std::string &a) { return a == "x" ? 0 : a == "y" ? 1 : 2; }

/// The material maps with the PML included, the mode-source update and its CPML read them
void fill_full_material_maps()
{
	const auto &S = subdomain;
	std::vector<int> background(S.cells_y + 1, sp.medium[0]);
	for (int j = 0; j <= S.cells_y && background_y.size(); ++j)
		background[j] = background_y[std::clamp(j + S.j_for_uni, 0, int(background_y.size()) - 1)];
	fdtd::materials::material_map().extend_into_pml(S.inner_org_x, S.inner_org_y, S.inner_org_z, S.inner_dim_x, S.inner_dim_y, S.inner_dim_z, S.cpml_only, background);
}
} // namespace

void fdtd::sources::setup_mode_source()
{
	fill_full_material_maps();

	const auto &m = mode_input;
	const std::string base = sp.input_dir;
	auto path = [&](const std::string &f) { return f.empty() || f[0] == '/' ? f : base + "/" + f; };
	MpbField e = read_mpb(path(m.efield), m.unit_length);
	MpbField h = read_mpb(path(m.hfield), m.unit_length);
	if (e.nx != h.nx || e.ny != h.ny)
		throw std::runtime_error("mode source: the E and H files have different grids");

	auto &P = plane;
	P.n = axis_index(m.propagation_axis);
	P.u = (P.n + 1) % 3, P.v = (P.n + 2) % 3;
	// FDTD axes of MPB's x and y
	int ax = m.mpb_x_axis.empty() ? (P.n == 0 ? 2 : 0) : axis_index(m.mpb_x_axis);
	int ay = m.mpb_y_axis.empty() ? (P.n == 1 ? 2 : 1) : axis_index(m.mpb_y_axis);

	std::vector<double> ex;
	std::vector<double> ey;
	std::vector<double> hx;
	std::vector<double> hy;
	make_real(e, ex, ey);
	make_real(h, hx, hy);
	// scale: |E_t| peaks at 1; H in A/m, MPB's H is in units of E: H_SI = H / eta0
	double top = 0;
	for (size_t n = 0; n < ex.size(); ++n)
		top = std::max(top, std::hypot(ex[n], ey[n]));
	if (!(top > 0))
		throw std::runtime_error("mode source: the E field of " + path(m.efield) + " is zero");
	const double se = 1 / top;
	const double sh = fdtd::constants::INVVACIMP / top;
	// MPB component along FDTD axis a, at MPB sample n
	auto mpb = [&](const std::vector<double> &cx, const std::vector<double> &cy, int a, size_t n) { return a == ax ? cx[n] : cy[n]; };
	// power along +n decides the sign of H, the files' propagation direction and axis handedness drop out
	double flow = 0;
	for (size_t n = 0; n < ex.size(); ++n)
		flow += mpb(ex, ey, P.u, n) * mpb(hx, hy, P.v, n) - mpb(ex, ey, P.v, n) * mpb(hx, hy, P.u, n);
	const double hs = flow < 0 ? -sh : sh;
	if (fdtd::grid::world_rank == ROOTPROCESS)
		fdtd::io::log("sources", "mode source: %s, %d x %d samples over %.4g x %.4g um, n_eff %.6g, along +%c at %.4g um; MPB x along %c, y along %c\n",
					  path(m.efield).c_str(), e.nx, e.ny, e.lx * 1e6, e.ly * 1e6, m.effective_index, "xyz"[P.n], m.position * 1e6, "xyz"[ax], "xyz"[ay]);

	if (fdtd::grid::world_rank == ROOTPROCESS)
	{
		// where the mode sits in MPB's cell: a quick check of the placement against the structure
		double w = 0;
		double cx = 0;
		double cy = 0;
		for (int i = 0; i < e.nx; ++i)
			for (int j = 0; j < e.ny; ++j)
			{
				const size_t n = size_t(i) * e.ny + j;
				const double a2 = ex[n] * ex[n] + ey[n] * ey[n];
				w += a2, cx += a2 * (-0.5 * e.lx + i * e.lx / e.nx), cy += a2 * (-0.5 * e.ly + j * e.ly / e.ny);
			}
		fdtd::io::log("sources", "mode source: |E_t|^2 centroid at MPB (x, y) = (%.4g, %.4g) um from the cell centre, i.e. %c %+.4g um, %c %+.4g um from the structure centre\n",
					  cx / w * 1e6, cy / w * 1e6, "xyz"[ax], (cx / w + std::array<double, 3>{m.offset_x, m.offset_y, m.offset_z}[ax]) * 1e6, "xyz"[ay],
					  (cy / w + std::array<double, 3>{m.offset_x, m.offset_y, m.offset_z}[ay]) * 1e6);
	}

	const auto &S = subdomain;
	const std::array<int, 3> start{S.start_x, S.start_y, S.start_z};
	const std::array<int, 3> cells{S.cells_x, S.cells_y, S.cells_z};
	const std::array<double, 3> d{lat.d[0], lat.d[1], lat.d[2]};
	const std::array<int, 3> centre{sp.centre[0] + lat.cpml[0], sp.centre[1] + lat.cpml[1],
									sp.centre[2] + lat.cpml[2]};
	const std::array<double, 3> offset{m.offset_x, m.offset_y, m.offset_z};
	const int p = centre[P.n] + int(std::lround(m.position / d[P.n]));
	P.e_local = (p - start[P.n] >= 0 && p - start[P.n] <= cells[P.n]) ? p - start[P.n] : -1;
	P.h_local = (p - 1 - start[P.n] >= 0 && p - 1 - start[P.n] <= cells[P.n]) ? p - 1 - start[P.n] : -1;
	P.active = P.e_local >= 0 || P.h_local >= 0;
	P.delay = 0.5 * d[P.n] * m.effective_index / fdtd::constants::C0;
	if (!P.active)
		return;

	P.nu = cells[P.u] + 1, P.nv = cells[P.v] + 1;
	const size_t N = size_t(P.nu) * P.nv;
	P.eu.assign(N, 0), P.ev.assign(N, 0), P.hu_c.assign(N, 0), P.hv_c.assign(N, 0), P.eu_c.assign(N, 0), P.ev_c.assign(N, 0);
	// coordinate [m] of local index a along axis A, shifted by s cells, from the MPB cell centre
	auto coord = [&](int A, int a, double s) { return (start[A] + a + s - centre[A]) * d[A] - offset[A]; };
	// profile of field (cx, cy) along FDTD axis a at a point with FDTD transverse coordinates (cu, cv)
	auto sample = [&](const MpbField &f, const std::vector<double> &cx, const std::vector<double> &cy, int a, double cu, double cv) {
		const double px = ax == P.u ? cu : cv;
		const double py = ay == P.u ? cu : cv;
		return f.at(a == ax ? cx : cy, px, py);
	};
	const double ch = sp.dt / (fdtd::constants::MU0 * d[P.n]);
	ctensor<int> *full[3] = {&matind_x_full, &matind_y_full, &matind_z_full};
	auto factor = [&](int A, int au, int av) {
		std::array<int, 3> ijk{};
		ijk[P.n] = std::max(P.e_local, 0), ijk[P.u] = au, ijk[P.v] = av;
		double c = 0;
		return e_factor((*full[A])[ijk[0]][ijk[1]][ijk[2]], c) ? c / d[P.n] : 0.0;
	};
	for (int au = 0; au < P.nu; ++au)
		for (int av = 0; av < P.nv; ++av)
		{
			const size_t q = size_t(au) * P.nv + av;
			// Yee positions in the plane: E_u (u+1/2, v), E_v (u, v+1/2), H_u (u, v+1/2), H_v (u+1/2, v)
			P.eu[q] = se * sample(e, ex, ey, P.u, coord(P.u, au, 0.5), coord(P.v, av, 0));
			P.ev[q] = se * sample(e, ex, ey, P.v, coord(P.u, au, 0), coord(P.v, av, 0.5));
			const double hu = hs * sample(h, hx, hy, P.u, coord(P.u, au, 0), coord(P.v, av, 0.5));
			const double hv = hs * sample(h, hx, hy, P.v, coord(P.u, au, 0.5), coord(P.v, av, 0));
			P.hv_c[q] = factor(P.u, au, av) * hv; // E_u correction
			P.hu_c[q] = factor(P.v, au, av) * hu; // E_v correction
			P.eu_c[q] = ch * P.eu[q];             // H_v correction
			P.ev_c[q] = ch * P.ev[q];             // H_u correction
		}
}

void fdtd::sources::mode_source_h(const Signal &signal, double t)
{
	const auto &P = plane;
	if (P.h_local < 0)
		return;
	// E_inc on the E plane at the time the H update used: half a step before the H the E update
	// uses, and half a cell downstream, measured: -46 dB backward against -28 dB with +1/2 step
	const double s = signal.value(t - 0.5 * sp.dt - P.delay);
	auto H = h_fields();
	for (int au = 0; au < P.nu; ++au)
		for (int av = 0; av < P.nv; ++av)
		{
			const size_t q = size_t(au) * P.nv + av;
			cell(*H[P.v], P.h_local, au, av) += s * P.eu_c[q];
			cell(*H[P.u], P.h_local, au, av) -= s * P.ev_c[q];
		}
}

void fdtd::sources::mode_source_e(const Signal &signal, double t)
{
	const auto &P = plane;
	if (P.e_local < 0)
		return;
	// H_inc half a cell upstream: it leads the E plane by the delay
	const double s = signal.value(t);
	auto E = e_fields();
	for (int au = 0; au < P.nu; ++au)
		for (int av = 0; av < P.nv; ++av)
		{
			const size_t q = size_t(au) * P.nv + av;
			cell(*E[P.u], P.e_local, au, av) += s * P.hv_c[q];
			cell(*E[P.v], P.e_local, au, av) -= s * P.hu_c[q];
		}
}
