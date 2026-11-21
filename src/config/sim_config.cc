/** \file   config/sim_config.cc
 *  \brief  Schema, loading, validation and serialization of SimConfig.
 *
 *  The schema is declared ONCE in visit_sim_config(). Adding a setting means
 *  adding a member to config/sim_config.h and one line here; loading, unknown-key
 *  detection, the resolved dump and the template all follow automatically.
 **/
#include "config/sim_config.h"
#include <stdexcept>
#include "geometry/shapes.h"
#include "sources/kinds.h"
#include "config/errors.h"
#include "core/quadrature.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include "materials/slots.h"

namespace fdtd::config
{
namespace
{

// The schema

const std::vector<std::string> PULSES = {
	"raised_cosine", "raised_cosine_derivative", "continuous_wave", "cw_pulse", "rectangular",
	"gaussian", "gaussian_alt", "constant", "sinc"};
/// [geometry] shape names and the solver's structure codes (antenna_type)
const std::vector<std::pair<std::string, int>> SHAPES = {
	{"none", fdtd::geometry::shape::none},
	{"sphere", fdtd::geometry::shape::sphere},
	{"dipole", fdtd::geometry::shape::dipole},
	{"bowtie", fdtd::geometry::shape::bowtie},
	{"cylinder", fdtd::geometry::shape::cylinder},
	{"cylinder_shell", fdtd::geometry::shape::cylinder_shell},
	{"two_spheres", fdtd::geometry::shape::two_spheres},
	{"l_antenna", fdtd::geometry::shape::l_antenna},
	{"hole", fdtd::geometry::shape::hole},
	{"lateral_dipole", fdtd::geometry::shape::lateral_dipole},
	{"round_butterfly", fdtd::geometry::shape::round_butterfly},
	{"round_dipole", fdtd::geometry::shape::round_dipole},
	{"one_side_round_dipole", fdtd::geometry::shape::one_side_round_dipole},
	{"round_tripole", fdtd::geometry::shape::round_tripole},
	{"round_hybrid", fdtd::geometry::shape::round_hybrid},
	{"x_shape", fdtd::geometry::shape::x_shape},
	{"sphere_dimer", fdtd::geometry::shape::sphere_dimer},
	{"sphere_trimer", fdtd::geometry::shape::sphere_trimer},
	{"sphere_quadrumer", fdtd::geometry::shape::sphere_quadrumer},
	{"ellipse", fdtd::geometry::shape::ellipse},
	{"sphere_shell", fdtd::geometry::shape::sphere_shell},
	{"triangle", fdtd::geometry::shape::triangle},
	{"rounded_triangle", fdtd::geometry::shape::rounded_triangle},
	{"rounded_monopole", fdtd::geometry::shape::rounded_monopole},
	{"u_shape", fdtd::geometry::shape::u_shape},
	{"cylinder_3d", fdtd::geometry::shape::cylinder_3d},
	{"h_shape", fdtd::geometry::shape::h_shape},
	{"v_shape", fdtd::geometry::shape::v_shape},
	{"l_shape", fdtd::geometry::shape::l_shape},
	{"two_films_spacer", fdtd::geometry::shape::two_films_spacer},
	{"rounded_dipole", fdtd::geometry::shape::rounded_dipole},
	{"rounded_dipole_offset_planes", fdtd::geometry::shape::rounded_dipole_offset_planes},
	{"anapole", fdtd::geometry::shape::anapole},
};

const double BUTTERFLY_ZERO = 45;    ///< [deg] the builder draws the bars along the diagonals at 0
const double ROUND_DIPOLE_ZERO = 90; ///< [deg] the builder turns the top arm from +x toward +z; the key counts from +z

// Every shape has a section of its own, [<shape>], with keys named for that shape, the anapole's is
// AnapoleSection. Each key is stored in the GeometrySection member the builder reads
// src/app/structures.cc, create_structures, so the builders are untouched.
using G = GeometrySection;
const std::vector<ShapeSection> SHAPE_SECTIONS = {
	{"sphere", {{"radius", &G::width, "[m]"}}},
	{"sphere_shell", {{"radius", &G::width, "[m] radius of the core; the shell of [geometry] shell_material is 13 cells thicker"}}},
	{"two_spheres", {{"radius", &G::width, "[m]"}, {"separation", &G::shift_y, "[m] centre to centre, along y"}}},
	{"sphere_dimer", {{"radius", &G::width, "[m]"}, {"spacing", &G::length, "[m] centre to centre, along z"}}},
	{"sphere_trimer", {{"radius", &G::width, "[m]"}, {"spacing", &G::length, "[m] centre to centre, along z"}}},
	{"sphere_quadrumer", {{"radius", &G::width, "[m]"}, {"spacing", &G::length, "[m] centre to centre, along z"}}},
	{"cylinder", {{"radius", &G::width, "[m] cylinder along z"}, {"length", &G::length, "[m] along z; 0 = through the whole domain"}}},
	{"cylinder_3d", {{"radius", &G::width, "[m] disc in the x-z plane"}, {"height", &G::length, "[m] along y"}}},
	{"cylinder_shell", {{"radius_core", &G::width, "[m] radius of the core; both cylinders are along y, on the same axis"},
						{"core_length", &G::length, "[m] of the core along y; 0 = through the whole domain"},
						{"radius_shell", &G::length_2, "[m] outer radius of the shell of [geometry] shell_material; the core is drawn over it"},
						{"shell_length", &G::arm_length, "[m] of the shell along y; 0 = through the whole domain; shorter than the core leaves the core ends bare"}}},
	{"ellipse", {{"radius_x", &G::width, "[m] semi-axis along x"}, {"radius_z", &G::length, "[m] semi-axis along z"}, {"thickness", &G::thickness, "[m] along y"},
				 {"rotation", &G::rotation, "[deg] positive turns counter-clockwise in the top view, x right and z up"}}},
	{"hole", {{"radius", &G::width, "[m] of the hole"}, {"thickness", &G::thickness, "[m] of the film that fills the x-z plane, along y"}}},
	{"dipole", {{"width", &G::width, "[m] arm width, along x"}, {"thickness", &G::thickness, "[m] along y"}, {"gap", &G::gap, "[m] between the arms, along z"},
				{"top_end", &G::length, "[m] end of the top arm along +z, measured from the centre of the gap"},
				{"bottom_end", &G::length_2, "[m] end of the bottom arm along -z, measured from the centre of the gap"}}},
	{"bowtie", {{"width", &G::width, "[m] base of each triangle, along x"}, {"thickness", &G::thickness, "[m] along y"},
				{"length", &G::length, "[m] base to base, along z"}, {"gap", &G::gap, "[m] between the tips"}}},
	{"triangle", {{"width", &G::width, "[m] base, along x"}, {"thickness", &G::thickness, "[m] along y"}, {"length", &G::length, "[m] base to apex, along z"}}},
	{"rounded_triangle", {{"width", &G::width, "[m] outer size along x"}, {"thickness", &G::thickness, "[m] along y"}, {"length", &G::length, "[m] outer size along z"},
						  {"corner_radius", &G::arm_length, "[m]"}}},
	{"rounded_monopole", {{"width", &G::width, "[m] along x"}, {"thickness", &G::thickness, "[m] along y"}, {"length", &G::length, "[m] along z"},
						  {"corner_radius", &G::arm_length, "[m] of the four vertical edges"}}},
	{"rounded_dipole", {{"width", &G::width, "[m] arm width, along x"}, {"thickness", &G::thickness, "[m] along y"}, {"length", &G::length, "[m] of each arm, along z"},
						{"gap", &G::gap, "[m] between the arms"}, {"corner_radius", &G::arm_length, "[m] of the vertical edges"}}},
	{"rounded_dipole_offset_planes", {{"width", &G::width, "[m] arm width, along x"}, {"thickness", &G::thickness, "[m] along y"}, {"length", &G::length, "[m] of each arm, along z"},
									  {"gap", &G::gap, "[m] between the arms"}, {"corner_radius", &G::arm_length, "[m] of the vertical edges"}}},
	{"round_dipole", {{"width", &G::width, "[m] arm width"}, {"thickness", &G::thickness, "[m] along y"}, {"gap", &G::gap, "[m] between the arms"},
					  {"top_length", &G::length_2, "[m] length of the top arm"}, {"bottom_length", &G::arm_width, "[m] length of the bottom arm"},
					  {"corner_radius", &G::corner_radius, "[m] of every edge"}, {"rotation", &G::rotation, "[deg] 0: arms along z; positive turns counter-clockwise in the top view, x right and z up", ROUND_DIPOLE_ZERO}}},
	{"one_side_round_dipole", {{"width", &G::width, "[m] arm width; the outer ends are half discs of this diameter"}, {"thickness", &G::thickness, "[m] along y"},
							   {"gap", &G::gap, "[m] between the arms"},
							   {"top_end", &G::length, "[m] centre of the top end disc, measured from the centre of the gap"},
							   {"bottom_end", &G::arm_length, "[m] centre of the bottom end disc, measured from the centre of the gap"},
							   {"rotation", &G::rotation, "[deg] 0: arms along z; positive turns counter-clockwise in the top view, x right and z up", ROUND_DIPOLE_ZERO}}},
	{"round_hybrid", {{"width", &G::width, "[m] arm width"}, {"thickness", &G::thickness, "[m] along y"}, {"gap", &G::gap, "[m] between the arms, filled with the gap material"},
					  {"arm_length", &G::length_2, "[m] length of each arm"}, {"corner_radius", &G::corner_radius, "[m] of every edge"},
					  {"rotation", &G::rotation, "[deg] 0: arms along z; positive turns counter-clockwise in the top view, x right and z up", ROUND_DIPOLE_ZERO}},
	 // the builder takes the arm length as (length_2 + arm_width + width - gap) / 2
	 [](GeometrySection& g) { g.length_2 = 2 * g.length_2 + g.gap - g.width; g.arm_width = 0; }},
	{"u_shape", {{"width", &G::width, "[m] width of the three bars"}, {"thickness", &G::thickness, "[m] along y"}, {"length", &G::length, "[m] outer size along z"},
				 {"arm_length", &G::length_2, "[m] of the two arms along x, beyond the base bar"}, {"corner_radius", &G::arm_length, "[m]"}}},
	{"h_shape", {{"width", &G::width, "[m] width of the two bars, along x"}, {"thickness", &G::thickness, "[m] along y"}, {"length", &G::length, "[m] of the two bars, along z"},
				 {"span", &G::length_2, "[m] outer size along x"}, {"bridge_width", &G::gap, "[m] size of the bridge along z"}, {"corner_radius", &G::arm_length, "[m]"}}},
	{"v_shape", {{"width", &G::width, "[m] arm width"}, {"thickness", &G::thickness, "[m] along y"}, {"length", &G::length, "[m] arm length"},
				 {"corner_radius", &G::arm_length, "[m]"}}},
	{"l_shape", {{"width", &G::width, "[m] arm width"}, {"thickness", &G::thickness, "[m] along y"}, {"length", &G::length, "[m] arm length"},
				 {"corner_radius", &G::arm_length, "[m]"}}},
	{"two_films_spacer", {{"thickness", &G::thickness, "[m] of the film, along y"}}},
	{"l_antenna", {{"width", &G::width, "[m] width of the two bars"}, {"thickness", &G::thickness, "[m] along y"},
				   {"length_x", &G::length, "[m] of the bar along x at the -z end, from the outer corner to the centre of its round end"},
				   {"length_z", &G::length_2, "[m] of the bar along z at the -x end, from the outer corner to the centre of its round end"}}},
	{"lateral_dipole", {{"width", &G::width, "[m] width of each bar, along x"}, {"thickness", &G::thickness, "[m] along y"},
						{"gap", &G::gap, "[m] between the two bars, along x"},
						{"left_length", &G::length, "[m] of the bar at -x; it reaches toward +z"},
						{"right_length", &G::length_2, "[m] of the bar at +x; it reaches toward -z"},
						{"overlap", &G::arm_length, "[m] length along z over which the two bars face each other; the facing corners are cut at 45 deg"}}},
	{"round_butterfly", {{"width", &G::width, "[m] bar width"}, {"thickness", &G::thickness, "[m] along y"}, {"gap", &G::gap, "[m] between the two halves; each half is one z bar and one x bar joined by a diagonal neck"},
						 {"length_z", &G::length_2, "[m] of each bar along z, measured from width / 2 off the centre"}, {"length_x", &G::arm_width, "[m] of each bar along x, measured the same way"},
						 {"corner_radius", &G::corner_radius, "[m] of every edge and of the two arcs"}, {"rotation", &G::rotation, "[deg] 0: bars along x and z; positive turns counter-clockwise in the top view, x right and z up", BUTTERFLY_ZERO}}},
	{"round_tripole", {{"width", &G::width, "[m] bar width"}, {"thickness", &G::thickness, "[m] along y"},
					   {"center_length", &G::length, "[m] of the middle bar; the three bars are in line"}, {"top_length", &G::length_2, "[m] of the top bar, along +z at rotation 0"}, {"bottom_length", &G::arm_width, "[m] of the bottom bar"},
					   {"top_gap", &G::arm_length, "[m] between the middle and the top bar"}, {"bottom_gap", &G::gap, "[m] between the middle and the bottom bar"},
					   {"corner_radius", &G::corner_radius, "[m] of every edge"}, {"rotation", &G::rotation, "[deg] 0: bars along z; positive turns counter-clockwise in the top view, x right and z up", ROUND_DIPOLE_ZERO}}},
	{"x_shape", {{"thickness", &G::thickness, "[m] along y"},
				 {"fixed_length", &G::length_2, "[m] of the fixed bar, along x"}, {"fixed_width", &G::width, "[m] of the fixed bar"},
				 {"rotated_length", &G::length, "[m] of the second bar, which crosses the fixed one at its centre"}, {"rotated_width", &G::gap, "[m] of the second bar"},
				 {"corner_radius", &G::corner_radius, "[m] of every edge"},
				 {"angle", &G::rotation, "[deg] from the fixed bar to the second bar, counter-clockwise in the top view, x right and z up; 90 is a cross"}}},
};

const double ZERO = 0;
double& shape_param(SimConfig& c, const char* shape, const char* key) { return c.shape_params[shape][key]; }
const double& shape_param(const SimConfig& c, const char* shape, const char* key)
{
	auto s = c.shape_params.find(shape);
	if (s == c.shape_params.end())
		return ZERO;
	auto k = s->second.find(key);
	return k == s->second.end() ? ZERO : k->second;
}

std::vector<std::string> shape_names()
{
	std::vector<std::string> n;
	for (const auto& s : SHAPES)
		n.push_back(s.first);
	return n;
}

const std::vector<std::string> SOURCES = {"plane_wave", "dipole", "gaussian_beam", "wave_superposition", "mode"};

/** Every scalar setting, exactly once. V must provide
 *    field(section, key, T& ref, doc, required)
 *    choice(section, key, std::string& ref, choices, doc) */
/// [monitor.<name>]: one box of cells recorded in time or frequency
template <class M, class V>
void visit_monitor(M& m, V& v, const std::string& sect)
{
	const char* s = sect.c_str();
	v.choice(s, "type", m.type, {"frequency", "time"}, "frequency: DFT of E at every [spectrum] point; time: E every `step` steps");
	v.field(s, "center_x", m.center_x, "box centre [cells, PML excluded]", true);
	v.field(s, "center_y", m.center_y, "", true);
	v.field(s, "center_z", m.center_z, "", true);
	v.field(s, "size_x", m.size_x, "box size [cells]", true);
	v.field(s, "size_y", m.size_y, "", true);
	v.field(s, "size_z", m.size_z, "", true);
	v.field(s, "step", m.step, "time monitors: record every N steps");
}

/// [dipole.<name>]: one point dipole of a dipole source
template <class D, class V>
void visit_dipole(D& d, V& v, const std::string& sect)
{
	const char* s = sect.c_str();
	v.field(s, "x", d.x, "position [cells, PML excluded]", true);
	v.field(s, "y", d.y, "", true);
	v.field(s, "z", d.z, "", true);
	v.field(s, "e_x", d.e_x, "amplitude of the injected Ex, times the pulse");
	v.field(s, "e_y", d.e_y, "");
	v.field(s, "e_z", d.e_z, "");
}

/// [block.<name>]: one box of a layered or periodic structure
template <class B, class V>
void visit_block(B& b, V& v, const std::string& sect)
{
	const char* s = sect.c_str();
	v.field(s, "material", b.material_name, "material name", true);
	v.field(s, "x_min", b.x_min, "[m] from the structure centre; +-1 m: past the domain");
	v.field(s, "x_max", b.x_max, "[m]");
	v.field(s, "y_min", b.y_min, "[m]");
	v.field(s, "y_max", b.y_max, "[m]");
	v.field(s, "z_min", b.z_min, "[m]");
	v.field(s, "z_max", b.z_max, "[m]");
	v.field(s, "count", b.count, "copies");
	v.field(s, "pitch_x", b.pitch_x, "[m] between copies");
	v.field(s, "pitch_y", b.pitch_y, "[m]");
	v.field(s, "pitch_z", b.pitch_z, "[m]");
}

/// [wave.<name>]: one plane wave of a wave_superposition source
template <class P, class V>
void visit_wave(P& w, V& v, const std::string& sect)
{
	const char* s = sect.c_str();
	v.field(s, "direction_x", w.direction_x, "unit propagation vector", true);
	v.field(s, "direction_y", w.direction_y, "", true);
	v.field(s, "direction_z", w.direction_z, "", true);
	v.field(s, "e_x", w.e_x, "E amplitude");
	v.field(s, "e_y", w.e_y, "");
	v.field(s, "e_z", w.e_z, "");
	v.field(s, "h_x", w.h_x, "H amplitude in vacuum units, scaled by n of the first medium");
	v.field(s, "h_y", w.h_y, "");
	v.field(s, "h_z", w.h_z, "");
}

template <class C, class V>
void visit_sim_config(C& c, V& v)
{
	v.field("run", "wallclock_limit_hours", c.run.wallclock_limit_hours, "stop before this wall time [h]; 0 = never");
	v.field("run", "wallclock_check_step", c.run.wallclock_check_step, "step at which the wall-clock limit is extrapolated");
	v.field("run", "trace_interval", c.run.trace_interval, "print fields at the animation point every N steps; 0 = off");
	v.field("run", "checkpoint_count", c.run.checkpoint_count, "probes/trackers/Q flushed this many times per run");

	v.field("paths", "input_dir", c.paths.input_dir, "data files: mode profiles, electron temperature; CLI -i");
	v.field("paths", "output_dir", c.paths.output_dir, "output root; CLI -o");

	v.field("grid", "nx", c.grid.nx, "interior cells along x, without the PML", true);
	v.field("grid", "ny", c.grid.ny, "interior cells along y", true);
	v.field("grid", "nz", c.grid.nz, "interior cells along z", true);
	v.field("grid", "dx", c.grid.dx, "cell size along x [m]", true);
	v.field("grid", "dy", c.grid.dy, "cell size along y [m]", true);
	v.field("grid", "dz", c.grid.dz, "cell size along z [m]", true);
	v.field("grid", "pml_x", c.grid.pml_x, "CPML cells per x face");
	v.field("grid", "pml_y", c.grid.pml_y, "CPML cells per y face");
	v.field("grid", "pml_z", c.grid.pml_z, "CPML cells per z face");
	v.field("grid", "time_steps", c.grid.time_steps, "FDTD time steps");
	v.field("grid", "courant", c.grid.courant, "Courant number");
	v.field("grid", "origin_x", c.grid.origin_x, "structure origin [cells]; -1 = centre");
	v.field("grid", "origin_y", c.grid.origin_y, "structure origin [cells]; -1 = centre");
	v.field("grid", "origin_z", c.grid.origin_z, "structure origin [cells]; -1 = centre");
	v.field("grid", "periodic_x", c.grid.periodic_x, "periodic boundary along x");
	v.field("grid", "periodic_z", c.grid.periodic_z, "periodic boundary along z");
	v.field("grid", "staircase", c.grid.staircase, "true: staircased material assignment; false: per Yee component");

	v.field("parallel", "ranks_x", c.parallel.ranks_x, "MPI ranks along x");
	v.field("parallel", "ranks_y", c.parallel.ranks_y, "MPI ranks along y");
	v.field("parallel", "ranks_z", c.parallel.ranks_z, "MPI ranks along z");
	v.field("parallel", "chunks_x", c.parallel.chunks_x, "chunk start indices along x, e.g. \"0 40 90\"; empty = equal split");
	v.field("parallel", "chunks_y", c.parallel.chunks_y, "chunk start indices along y");
	v.field("parallel", "chunks_z", c.parallel.chunks_z, "chunk start indices along z");

	v.field("boxes", "tfsf_x", c.boxes.tfsf_x, "TF/SF box margin from the PML [cells]");
	v.field("boxes", "tfsf_y", c.boxes.tfsf_y, "");
	v.field("boxes", "tfsf_z", c.boxes.tfsf_z, "");
	v.field("boxes", "scattering_x", c.boxes.scattering_x, "scattering-flux box margin [cells]; -1 = tfsf - 5");
	v.field("boxes", "scattering_y", c.boxes.scattering_y, "");
	v.field("boxes", "scattering_z", c.boxes.scattering_z, "");
	v.field("boxes", "absorption_x", c.boxes.absorption_x, "absorption-flux box margin [cells]; -1 = tfsf + 5");
	v.field("boxes", "absorption_y", c.boxes.absorption_y, "");
	v.field("boxes", "absorption_z", c.boxes.absorption_z, "");

	v.choice("source", "type", c.source.type, SOURCES, "source; its keys are in the section of the same name");
	v.choice("source", "signal", c.source.signal, PULSES, "time signal; its keys are in [signal.<name>]");
	v.field("source", "peak_field", c.source.peak_field, "[V/m]");
	v.field("source", "allow_truncated_signal", c.source.allow_truncated_signal, "run even when the time steps do not cover the signal");
	v.field("source", "aux_steps", c.source.aux_steps, "1D auxiliary grid of plane_wave and wave_superposition: steps; -1 = [grid] time_steps");
	v.field("source", "aux_length_factor", c.source.aux_length_factor, "1D auxiliary grid: length as a multiple of the 3D y extent");

	v.choice("plane_wave", "propagation_axis", c.plane_wave.propagation_axis, {"y"}, "the wave travels along + this axis");
	v.field("plane_wave", "propagation_angle", c.plane_wave.propagation_angle, "[deg] only 0; oblique incidence is type = wave_superposition");
	v.choice("plane_wave", "polarization_axis", c.plane_wave.polarization_axis, {"x", "z", "rhcp", "lhcp"}, "axis of E, or circular polarisation");
	v.field("plane_wave", "polarization_angle", c.plane_wave.polarization_angle, "[deg] turns E from polarization_axis toward the other transverse axis");

	v.field("wave_superposition", "derivative", c.wave_superposition.derivative, "inject the time derivative of the signal");
	v.field("wave_superposition", "dispersion_correction", c.wave_superposition.dispersion_correction, "correct each wave for grid dispersion");
	v.field("wave_superposition", "layered_simulation", c.wave_superposition.layered_simulation, "run the 1D simulations of the waves even when the three layers are one medium");

	v.field("signal.raised_cosine", "max_frequency", c.signal.raised_cosine.max_frequency, "bandwidth [Hz]");
	v.field("signal.raised_cosine_derivative", "max_frequency", c.signal.raised_cosine_derivative.max_frequency, "bandwidth [Hz]");
	v.field("signal.rectangular", "max_frequency", c.signal.rectangular.max_frequency, "bandwidth [Hz]");
	v.field("signal.gaussian", "center_wavelength", c.signal.gaussian.center_wavelength, "carrier [m]");
	v.field("signal.gaussian", "width", c.signal.gaussian.width, "[s]");
	v.field("signal.gaussian", "chirp", c.signal.gaussian.chirp, "linear chirp");
	v.field("signal.gaussian", "chirp_quadratic", c.signal.gaussian.chirp_quadratic, "quadratic chirp");
	v.field("signal.gaussian_alt", "center_wavelength", c.signal.gaussian_alt.center_wavelength, "carrier [m]");
	v.field("signal.gaussian_alt", "width", c.signal.gaussian_alt.width, "[s]");
	v.field("mode", "efield", c.mode.efield, "MPB E-field file of the band, e.k01.b01.h5, relative to [paths] input_dir");
	v.field("mode", "hfield", c.mode.hfield, "MPB H-field file of the band, h.k01.b01.h5");
	v.field("mode", "unit_length", c.mode.unit_length, "[m] MPB's length unit a");
	v.field("mode", "effective_index", c.mode.effective_index, "n_eff = k / frequency of the band in MPB units");
	v.choice("mode", "propagation_axis", c.mode.propagation_axis, {"x", "y", "z"}, "the mode travels along + this axis");
	v.field("mode", "position", c.mode.position, "[m] injection plane along the propagation axis, from the structure centre");
	v.field("mode", "mpb_x_axis", c.mode.mpb_x_axis, "FDTD axis of MPB's x; empty: z for propagation along x, x otherwise");
	v.field("mode", "mpb_y_axis", c.mode.mpb_y_axis, "FDTD axis of MPB's y; empty: y for propagation along x or z, z for y");
	v.field("mode", "offset_x", c.mode.offset_x, "[m] MPB cell centre from the structure centre");
	v.field("mode", "offset_y", c.mode.offset_y, "[m]");
	v.field("mode", "offset_z", c.mode.offset_z, "[m]");

	v.field("signal.continuous_wave", "settle_periods", c.signal.continuous_wave.settle_periods, "the DFT starts after this many periods of the [spectrum] centre frequency");
	v.field("signal.cw_pulse", "periods", c.signal.cw_pulse.periods, "signal length [periods of the [spectrum] centre frequency]");

	v.field("gaussian_beam", "theta_steps", c.gaussian_beam.theta_steps, "");
	v.field("gaussian_beam", "phi_steps", c.gaussian_beam.phi_steps, "");
	v.field("gaussian_beam", "theta_max", c.gaussian_beam.theta_max, "[deg]");
	v.field("gaussian_beam", "phi_max", c.gaussian_beam.phi_max, "[deg]");
	v.field("gaussian_beam", "waist_factor", c.gaussian_beam.waist_factor, "beam waist = waist_factor * sin(theta_max) * focal length");
	v.field("gaussian_beam", "focal_factor", c.gaussian_beam.focal_factor, "focal length = focal_factor * distance from the grid origin to the structure centre");

	v.field("spectrum", "lambda_min", c.spectrum.lambda_min, "shortest DFT wavelength [m]", true);
	v.field("spectrum", "lambda_max", c.spectrum.lambda_max, "longest DFT wavelength [m]", true);
	v.field("spectrum", "points", c.spectrum.points, "samples in [lambda_min, lambda_max]; 1 needs lambda_min = lambda_max");
	v.choice("spectrum", "spacing", c.spectrum.spacing, {"wavelength", "frequency", "gauss"},
			 "uniform step in wavelength or frequency, or gauss: Gauss-Legendre nodes in frequency, the most accurate band integrals for the points");
	v.choice("spectrum", "harmonic_band", c.spectrum.harmonic_band, {"none", "second", "third"}, "extra DFT band at 2w / 3w of every point");

	v.field("output", "dft_xy", c.output.dft_xy, "frequency-domain E on the XY plot plane: fields/E_<wl>nm_XY.h5");
	v.field("output", "dft_xz", c.output.dft_xz, "frequency-domain E on the XZ plot plane");
	v.field("output", "dft_yz", c.output.dft_yz, "frequency-domain E on the YZ plot plane");
	v.field("output", "dft_magnetic", c.output.dft_magnetic, "also H on the enabled planes: fields/H_<wl>nm_<plane>.h5");
	v.field("output", "volume_dft", c.output.volume_dft, "full-volume E/H DFT");
	v.field("output", "frames_xy", c.output.frames_xy, "time-domain XY frame period [steps]; 0 = off");
	v.field("output", "frames_xz", c.output.frames_xz, "");
	v.field("output", "frames_yz", c.output.frames_yz, "");
	v.field("output", "frames_3d", c.output.frames_3d, "");
	v.field("output", "phase_normalize", c.output.phase_normalize, "normalize the DFT phase to the source");
	v.field("output", "xdmf", c.output.xdmf, "write a .xmf next to every HDF5 output, for ParaView / VisIt");
	v.field("output", "sigma_absorption", c.output.sigma_absorption, "absorption from the sigma|E|^2 volume integral");

	v.field("multipoles", "exact", c.multipoles.exact, "exact multipole decomposition of Alaee et al. 2018");
	v.field("multipoles", "lwa", c.multipoles.lwa, "long-wavelength multipoles, the kr -> 0 limit of the same integrals; any material model");
	v.field("multipoles", "sphere_radius", c.multipoles.sphere_radius, "legacy sphere-sampled multipoles: radius [cells]; 0 = off");
	v.field("multipoles", "sphere_theta", c.multipoles.sphere_theta, "");
	v.field("multipoles", "sphere_phi", c.multipoles.sphere_phi, "");

	v.choice("geometry", "shape", c.geometry.shape, shape_names(), "structure");
	v.field("geometry", "material", c.geometry.material, "material of the structure, by name: [material.<name>] or a catalog file");
	v.field("geometry", "first_medium", c.geometry.first_medium, "entry layer material");
	v.field("geometry", "middle_medium", c.geometry.middle_medium, "middle layer material");
	v.field("geometry", "last_medium", c.geometry.last_medium, "exit layer material");
	v.field("geometry", "middle_thickness", c.geometry.middle_thickness, "[m]");
	v.field("geometry", "shift_y", c.geometry.shift_y, "structure offset along y [m]");
	v.field("geometry", "shell_material", c.geometry.shell_material, "material of the shell, sphere_shell and cylinder_shell");

	for (const auto& sh : SHAPE_SECTIONS)
		for (const auto& k : sh.keys)
			v.field(sh.name, k.key, shape_param(c, sh.name, k.key), k.doc);

	v.field("materials", "catalog", c.materials.catalog, "directories of material files, ':' between them, relative to this file; then $FDTD_MATERIALS, then materials/ next to the executable");

	v.field("anapole", "patch_edge", c.anapole.patch_edge, "[m] edge of the square patch");
	v.field("anapole", "thickness", c.anapole.thickness, "[m] film thickness along y");
	v.field("anapole", "corner_radius", c.anapole.corner_radius, "[m] fillet radius of the patch corners");
	v.field("anapole", "lobe_length", c.anapole.lobe_length, "length of the lobe curve, in units of patch_edge");
	v.field("anapole", "lobe_height", c.anapole.lobe_height, "height coefficient of the lobe curve, in units of patch_edge");
	v.field("anapole", "lobe_exponent", c.anapole.lobe_exponent, "exponent of the lobe curve");
	v.field("anapole", "lobe_offset", c.anapole.lobe_offset, "where the lobe curve starts, from the patch centre, in units of patch_edge");
	v.field("anapole", "gap_length", c.anapole.gap_length, "[m] gap between the lobes, along the lobe axis");
	v.field("anapole", "gap_width", c.anapole.gap_width, "[m] gap across the lobe axis");
	v.field("anapole", "rotation", c.anapole.rotation, "[deg] in-plane rotation, x-z");
	v.field("anapole", "aperture_material", c.anapole_fill.aperture, "fill of the whole aperture, lobes and gap, a material name; unset = background");
	v.field("anapole", "gap_material", c.anapole_fill.gap, "fill of the gap alone, a material name; unset = as the rest of the aperture");
	v.field("anapole", "lobe_points", c.anapole.lobe_points, "polygon vertices per lobe edge");

	v.field("nonlinear", "shg", c.nonlinear.shg, "second-harmonic generation in every material with a chi2 tensor");
	v.field("nonlinear", "thg", c.nonlinear.thg, "third-harmonic generation in every material with a chi3 tensor");
	v.choice("nonlinear", "update", c.nonlinear.update, {"field", "polarization"}, "how the nonlinearity enters the update");
	v.field("nonlinear", "graphene_optimization", c.nonlinear.graphene_optimization, "");

	v.choice("electrons", "model", c.electrons.model, {"local", "nonlocal", "gnor", "hydrodynamic"}, "free-electron response; the hydrodynamic keys are in [hydrodynamic]");
	v.field("electrons", "fermi_velocity", c.electrons.fermi_velocity, "factor on the fermi_velocity of every material, nonlocal, gnor, hydrodynamic");

	v.field("hydrodynamic", "pressure", c.hydrodynamic.pressure, "quantum pressure term");
	v.field("hydrodynamic", "convection", c.hydrodynamic.convection, "convective term");
	v.field("hydrodynamic", "magnetic", c.hydrodynamic.magnetic, "magnetic Lorentz-force term");
	v.field("hydrodynamic", "split_nonlinear", c.hydrodynamic.split_nonlinear, "separate linear / nonlinear DFTs");

	v.choice("temperature", "model", c.temperature.model, {"none", "two_temperature", "load", "pump", "probe"},
			 "electron temperature: two-temperature model, load <input_dir>/electron_temperature.h5, pump-probe");

	v.field("two_temperature", "save_frames", c.two_temperature.save_frames, "write temperature frames of model = two_temperature");
	v.field("two_temperature", "rk_stages", c.two_temperature.rk_stages, "Runge-Kutta stages");

	v.field("pump_probe", "delay", c.pump_probe.delay, "[s]");
}

// Value codecs

void parse_value(const io::Ini& ini, const std::string& s, const std::string& k, int& out)
{
	out = io::Ini::to_integer(s, k, ini.get(s, k));
}
void parse_value(const io::Ini& ini, const std::string& s, const std::string& k, double& out)
{
	out = io::Ini::to_double(s, k, ini.get(s, k));
}
void parse_value(const io::Ini& ini, const std::string& s, const std::string& k, bool& out)
{
	out = ini.get_bool(s, k);
}
void parse_value(const io::Ini& ini, const std::string& s, const std::string& k, std::string& out)
{
	out = io::Ini::unquote(ini.get(s, k));
}

std::string format_value(int v) { return std::to_string(v); }
std::string format_value(bool v) { return v ? "true" : "false"; }
std::string format_value(const std::string& v) { return v; }
std::string format_value(double v)
{
	char buf[64];
	auto res = std::to_chars(buf, buf + sizeof(buf), v); // shortest round-trip representation
	return std::string(buf, res.ptr);
}

std::string lower(std::string s) { return io::Ini::lower(std::move(s)); }

// Visitors

struct Loader
{
	const io::Ini& ini;
	std::set<std::pair<std::string, std::string>> known;
	std::vector<std::string> missing;

	template <class T>
	void field(const char* s, const char* k, T& ref, const char*, bool required = false)
	{
		known.insert({s, k});
		if (ini.is_key(s, k))
			parse_value(ini, s, k, ref);
		else if (required)
			missing.push_back(std::string("[") + s + "] " + k);
	}

	void choice(const char* s, const char* k, std::string& ref, const std::vector<std::string>& choices, const char* doc)
	{
		field(s, k, ref, doc);
		ref = lower(ref);
		check_choice(s, k, ref, choices);
	}

	static void check_choice(const char* s, const char* k, const std::string& ref, const std::vector<std::string>& choices)
	{
		if (std::find(choices.begin(), choices.end(), ref) == choices.end())
		{
			std::string all;
			for (const auto& c : choices)
				all += (all.empty() ? "" : ", ") + c;
			throw syntax_error(s, k, "invalid value '" + ref + "', expected one of: " + all);
		}
	}
};

struct Writer
{
	std::ostringstream out;
	std::string current;
	bool with_docs = false;

	void header(const std::string& s)
	{
		if (s == current)
			return;
		out << (current.empty() ? "" : "\n") << "[" << s << "]\n";
		current = s;
	}

	template <class T>
	void field(const char* s, const char* k, const T& ref, const char* doc, bool required = false)
	{
		header(s);
		if (with_docs && (required || (doc && *doc)))
			out << "; " << (required ? "required: " : "") << (doc ? doc : "") << "\n";
		out << k << " = " << format_value(ref) << "\n";
	}

	void choice(const char* s, const char* k, const std::string& ref, const std::vector<std::string>& choices, const char* doc)
	{
		header(s);
		if (with_docs)
		{
			out << "; " << (doc && *doc ? std::string(doc) + ": " : "");
			for (size_t i = 0; i < choices.size(); ++i)
				out << (i ? " | " : "") << choices[i];
			out << "\n";
		}
		out << k << " = " << ref << "\n";
	}

};

} // namespace

// Public api

namespace
{
int lookup(const std::map<std::string, int>& table, const std::string& name, const char* what)
{
	auto it = table.find(lower(name));
	if (it == table.end())
		throw validation_error(std::string("unknown ") + what + " '" + name + "'");
	return it->second;
}
} // namespace

int pulse_code(const std::string& name)
{
	return lookup({{"raised_cosine", fdtd::sources::pulse::raised_cosine}, {"raised_cosine_derivative", fdtd::sources::pulse::raised_cosine_derivative},
				   {"continuous_wave", fdtd::sources::pulse::continuous_wave}, {"cw_pulse", fdtd::sources::pulse::cw_pulse}, {"rectangular", fdtd::sources::pulse::rectangular},
				   {"gaussian", fdtd::sources::pulse::gaussian}, {"gaussian_alt", fdtd::sources::pulse::gaussian_alt}, {"constant", fdtd::sources::pulse::constant},
				   {"sinc", fdtd::sources::pulse::sinc}},
				  name, "[source] signal");
}

std::vector<int> tfsf_planes(const SimConfig& c, int axis)
{
	const int n[3] = {c.grid.nx, c.grid.ny, c.grid.nz};
	const int pml[3] = {c.grid.pml_x, c.grid.pml_y, c.grid.pml_z};
	const int tfsf[3] = {c.boxes.tfsf_x, c.boxes.tfsf_y, c.boxes.tfsf_z};
	const bool periodic[3] = {c.grid.periodic_x, false, c.grid.periodic_z};
	if (periodic[axis])
		return {};
	return {pml[axis] + tfsf[axis], pml[axis] + n[axis] - tfsf[axis]};
}

std::vector<int> rank_starts(const char* axis, int total, int ranks, const std::string& chunks,
							const std::vector<int>& forbidden, std::string* note)
{
	const std::string key = std::string("chunks_") + axis;
	auto on_plane = [&](int b) { return std::find(forbidden.begin(), forbidden.end(), b) != forbidden.end(); };
	if (ranks < 1 || 2 * ranks > total)
		throw validation_error("[parallel] ranks_" + std::string(axis) + " = " + std::to_string(ranks) +
									": every rank needs at least 2 of the " + std::to_string(total) +
									" cells along " + axis + ", PML included; use at most " + std::to_string(total / 2));
	std::vector<int> start;
	if (!chunks.empty())
	{
		start = parse_int_list("parallel", key, chunks);
		for (size_t r = 1; r < start.size(); ++r)
			if (on_plane(start[r]))
				throw validation_error("[parallel] " + key + ": the rank boundary at cell " + std::to_string(start[r]) +
											" lies on the TF/SF plane along " + axis +
											"; move it by one cell, to " + std::to_string(start[r] + 1));
		return start;
	}
	const int base = total / ranks;
	const int extra = total % ranks;
	start.assign(ranks, 0);
	for (int r = 1; r < ranks; ++r)
		start[r] = start[r - 1] + base + (r - 1 < extra ? 1 : 0);
	auto size_ok = [&](int r) {
		const int end = r + 1 < ranks ? start[r + 1] : total;
		return end - start[r] >= 2;
	};
	for (int r = 1; r < ranks; ++r)
	{
		if (!on_plane(start[r]))
			continue;
		const int was = start[r];
		bool moved = false;
		for (int shift : {+1, -1})
		{
			start[r] = was + shift;
			if (!on_plane(start[r]) && size_ok(r - 1) && size_ok(r))
			{
				moved = true;
				break;
			}
		}
		if (!moved)
			throw validation_error(std::string("[parallel] no valid split of ") + std::to_string(total) + " cells into " +
										std::to_string(ranks) + " ranks along " + axis +
										" keeps the rank boundaries off the TF/SF planes; change ranks_" + axis);
		if (note)
			*note += std::string(*note == "" ? "" : "\n") + "[parallel] the " + axis + " rank boundary at cell " +
					 std::to_string(was) + " fell on the TF/SF plane; moved to " + std::to_string(start[r]);
	}
	return start;
}

int shape_code(const std::string& name)
{
	for (const auto& s : SHAPES)
		if (s.first == lower(name))
			return s.second;
	throw validation_error("unknown [geometry] shape '" + name + "'");
}

int source_code(const std::string& name)
{
	return lookup({{"plane_wave", fdtd::sources::kind::plane_wave}, {"dipole", fdtd::sources::kind::dipole}, {"gaussian_beam", fdtd::sources::kind::gaussian_beam},
				   {"wave_superposition", fdtd::sources::kind::wave_superposition}, {"mode", fdtd::sources::kind::mode}},
				  name, "[source] type");
}

int dft_plane_mask_bits(const OutputSection& o)
{
	const int planes = (o.dft_yz ? 1 : 0) | (o.dft_xz ? 2 : 0) | (o.dft_xy ? 4 : 0);
	return planes ? planes | (o.dft_magnetic ? 8 : 0) : 0;
}

int harmonic_order(const std::string& name)
{
	return lookup({{"none", 1}, {"second", 2}, {"third", 3}}, name, "[spectrum] harmonic_band");
}

int process_harmonic(const NonlinearSection& n)
{
	return n.shg ? 2 : n.thg ? 3 : 1;
}

std::vector<double> spectrum_weights(const SpectrumSection& s, double c0, double two_pi)
{
	const int n = s.points;
	std::vector<double> q(n, 0.0);
	if (n < 2)
		return q;
	if (lower(s.spacing) == "gauss")
	{
		std::vector<double> x;
		fdtd::core::gauss_legendre(n, x, q);
		const double half = 0.5 * (c0 / s.lambda_min - c0 / s.lambda_max);
		for (double& v : q)
			v *= half;
		return q;
	}
	const std::vector<double> w = spectrum_omegas(s, c0, two_pi);
	for (int m = 0; m + 1 < n; ++m)
	{
		const double half = 0.5 * (w[m + 1] - w[m]) / two_pi;
		q[m] += half;
		q[m + 1] += half;
	}
	return q;
}

std::vector<double> spectrum_omegas(const SpectrumSection& s, double c0, double two_pi)
{
	std::vector<double> w;
	const int n = s.points;
	std::vector<double> nodes;
	std::vector<double> unused;
	if (n > 1 && lower(s.spacing) == "gauss")
		fdtd::core::gauss_legendre(n, nodes, unused);
	for (int m = 0; m < n; ++m)
	{
		if (n == 1)
			w.push_back(two_pi * (c0 / s.lambda_min));
		else if (lower(s.spacing) == "gauss")
		{
			const double f0 = c0 / s.lambda_max;
			const double f1 = c0 / s.lambda_min;
			w.push_back(two_pi * (0.5 * (f0 + f1) + 0.5 * (f1 - f0) * nodes[m]));
		}
		else if (lower(s.spacing) == "frequency")
		{
			const double f0 = c0 / s.lambda_max;
			const double f1 = c0 / s.lambda_min;
			w.push_back(two_pi * (f0 + m * (f1 - f0) / (n - 1)));
		}
		else // uniform in wavelength, ascending frequency
			w.push_back(two_pi * c0 / (s.lambda_max - m * (s.lambda_max - s.lambda_min) / (n - 1)));
	}
	const int k = harmonic_order(s.harmonic_band);
	if (k > 1)
		for (int m = 0; m < n; ++m)
			w.push_back(k * w[m]);
	return w;
}

std::vector<int> parse_int_list(const std::string& sect, const std::string& key, const std::string& raw)
{
	std::vector<int> out;
	std::istringstream in(raw);
	std::string tok;
	while (in >> tok)
		out.push_back(io::Ini::to_integer(sect, key, tok));
	return out;
}

bool split_override(const std::string& spec, std::string& section, std::string& key, std::string& value)
{
	const auto eq = spec.find('=');
	if (eq == std::string::npos)
		return false;
	std::string path = spec.substr(0, eq);
	value = spec.substr(eq + 1);
	io::Ini::trim(path);
	io::Ini::trim(value);
	// split at the LAST dot so that dynamic sections, "probe.top.x", work
	const auto dot = path.rfind('.');
	if (dot == std::string::npos || dot == 0 || dot + 1 == path.size())
		return false;
	section = path.substr(0, dot);
	key = path.substr(dot + 1);
	return true;
}

io::Ini build_config(const std::string& ini_text, const std::string& name, const std::vector<std::string>& overrides,
						const std::string& input_dir, const std::string& output_dir)
{
	io::Ini ini;
	ini.load_string(ini_text, name.empty() ? "<config>" : name);
	for (const auto& spec : overrides)
	{
		std::string s;
		std::string k;
		std::string v;
		if (!split_override(spec, s, k, v))
			throw std::invalid_argument("override expects SECTION.KEY=VALUE, got '" + spec + "'");
		ini.override(s, {{k, v}});
	}
	if (!input_dir.empty())
		ini.override("paths", {{"input_dir", input_dir}});
	if (!output_dir.empty())
		ini.override("paths", {{"output_dir", output_dir}});
	return ini;
}

SimConfig load_sim_config(const io::Ini& ini)
{
	SimConfig c;
	Loader loader{ini, {}, {}};
	visit_sim_config(c, loader);

	for (const auto& sect : ini.sections())
	{
		if (sect.rfind("probe.", 0) != 0)
			continue;
		ProbePoint p;
		p.name = sect.substr(6);
		for (const char* k : {"x", "y", "z"})
			if (!ini.is_key(sect, k))
				loader.missing.push_back("[" + sect + "] " + k);
		p.x = ini.get_int(sect, "x");
		p.y = ini.get_int(sect, "y");
		p.z = ini.get_int(sect, "z");
		c.probes.push_back(p);
	}
	for (const auto& sect : ini.sections())
	{
		if (sect.rfind("monitor.", 0) == 0)
		{
			MonitorBox m;
			m.name = sect.substr(8);
			visit_monitor(m, loader, sect);
			c.monitors.push_back(m);
		}
		else if (sect.rfind("dipole.", 0) == 0)
		{
			DipoleSpec d;
			d.name = sect.substr(7);
			visit_dipole(d, loader, sect);
			c.dipoles.push_back(d);
		}
		else if (sect.rfind("wave.", 0) == 0)
		{
			WaveSpec w;
			w.name = sect.substr(5);
			visit_wave(w, loader, sect);
			c.waves.push_back(w);
		}
	}

	for (const auto& sect : ini.sections_in_file_order()) // drawn in file order: a later block overwrites
		if (sect.rfind("block.", 0) == 0)
		{
			BlockSpec b;
			b.name = sect.substr(6);
			visit_block(b, loader, sect);
			c.blocks.push_back(b);
		}

	// unknown keys are errors: a typo must never be silently ignored
	std::vector<std::string> unknown;
	for (const auto& sect : ini.sections())
	{
		const bool is_probe = sect.rfind("probe.", 0) == 0;
		if (sect.rfind("chi2.", 0) == 0 || sect.rfind("chi3.", 0) == 0 || sect.rfind("material.", 0) == 0)
			continue; // the materials check their own keys, materials::definition
		for (const auto& key : ini.keys(sect))
		{
			if (is_probe && (key == "x" || key == "y" || key == "z"))
				continue;
			if (!loader.known.count({sect, key}))
				unknown.push_back("[" + sect + "] " + key);
		}
	}
	if (!unknown.empty() || !loader.missing.empty())
	{
		std::string msg;
		for (const auto& u : unknown)
			msg += "\n  unknown key:  " + u;
		for (const auto& m : loader.missing)
			msg += "\n  missing key:  " + m;
		throw validation_error("configuration errors:" + msg);
	}

	// the materials by name, in slot order; every material section of the file is checked
	for (std::string* name : {&c.geometry.material, &c.geometry.first_medium, &c.geometry.middle_medium, &c.geometry.last_medium,
							  &c.geometry.shell_material, &c.anapole_fill.aperture, &c.anapole_fill.gap})
		*name = lower(*name);
	for (auto& b : c.blocks)
		b.material_name = lower(b.material_name);
	try
	{
		for (const auto& sect : ini.sections())
			if (sect.rfind("material.", 0) == 0 && sect.find('.', 9) == std::string::npos)
				materials::definition(ini, sect.substr(9));
		for (const auto& name : material_names(c))
			c.definitions.push_back(materials::definition(ini, name));
	}
	catch (const std::invalid_argument& e)
	{
		throw validation_error(e.what());
	}

	// resolve derived defaults so the resolved dump shows real values
	auto def = [](int& v, int fallback) { if (v < 0) v = fallback; };
	def(c.boxes.scattering_x, c.boxes.tfsf_x - 5);
	def(c.boxes.scattering_y, c.boxes.tfsf_y - 5);
	def(c.boxes.scattering_z, c.boxes.tfsf_z - 5);
	def(c.boxes.absorption_x, c.boxes.tfsf_x + 5);
	def(c.boxes.absorption_y, c.boxes.tfsf_y + 5);
	def(c.boxes.absorption_z, c.boxes.tfsf_z + 5);
	def(c.grid.origin_x, c.grid.nx / 2);
	def(c.grid.origin_y, c.grid.ny / 2);
	def(c.grid.origin_z, c.grid.nz / 2);
	def(c.source.aux_steps, c.grid.time_steps);

	return c;
}

void validate_sim_config(const SimConfig& c)
{
	std::vector<std::string> err;
	auto require = [&err](bool ok, const std::string& msg) { if (!ok) err.push_back(msg); };

	const auto& g = c.grid;
	require(g.nx > 0 && g.ny > 0 && g.nz > 0, "[grid] nx, ny, nz must be > 0");
	require(g.dx > 0 && g.dy > 0 && g.dz > 0, "[grid] dx, dy, dz must be > 0");
	require(g.pml_x >= 0 && g.pml_y >= 0 && g.pml_z >= 0, "[grid] pml_* must be >= 0");
	require(g.time_steps > 0, "[grid] time_steps must be > 0");
	require(g.courant > 0 && g.courant <= 1, "[grid] courant must be in (0, 1]");
	require(c.parallel.ranks_x > 0 && c.parallel.ranks_y > 0 && c.parallel.ranks_z > 0, "[parallel] ranks_* must be >= 1");
	require(c.run.checkpoint_count >= 1, "[run] checkpoint_count must be >= 1");
	require(c.run.wallclock_limit_hours >= 0, "[run] wallclock_limit_hours must be >= 0");
	require(c.boxes.scattering_x >= 0 && c.boxes.scattering_y >= 0 && c.boxes.scattering_z >= 0,
			"[boxes] scattering_* must be >= 0: tfsf_* may be too small");

	// spectrum: each key means one thing
	const auto& s = c.spectrum;
	require(s.lambda_min > 0 && s.lambda_max > 0, "[spectrum] lambda_min and lambda_max must be > 0");
	require(s.points >= 1, "[spectrum] points must be >= 1");
	if (s.points == 1)
		require(s.lambda_min == s.lambda_max, "[spectrum] a single point needs lambda_min = lambda_max");
	else
		require(s.lambda_min < s.lambda_max, "[spectrum] lambda_min must be < lambda_max");

	const int band = harmonic_order(s.harmonic_band);
	const int process = process_harmonic(c.nonlinear);
	require(!(c.nonlinear.shg && c.nonlinear.thg), "[nonlinear] shg and thg cannot both be on");
	if (band > 1 && lower(c.electrons.model) == "local" && lower(c.temperature.model) == "none")
		require(process == band, "[spectrum] harmonic_band = " + s.harmonic_band + " requires [nonlinear] " + (band == 2 ? "shg" : "thg") + " = true");
	if (process > 1)
		require(band == process, std::string("[nonlinear] ") + (process == 2 ? "shg" : "thg") +
									 " computes the harmonic: set [spectrum] harmonic_band = " + (process == 2 ? "second" : "third"));
	{
		const std::string e = lower(c.electrons.model);
		const std::string t = lower(c.temperature.model);
		auto unused = [&](bool untouched, const char* section, const std::string& selector) {
			require(untouched, std::string("[") + section + "] has keys set, but " + selector + " does not use that section");
		};
		if (e == "local")
			require(c.electrons.fermi_velocity == ElectronsSection().fermi_velocity, "[electrons] fermi_velocity needs model = nonlocal, gnor or hydrodynamic");
		if (e != "hydrodynamic")
			unused(c.hydrodynamic == HydrodynamicSection(), "hydrodynamic", "[electrons] model = " + e);
		const bool ttm = t == "two_temperature" || t == "pump" || t == "probe";
		if (!ttm)
			unused(c.two_temperature == TwoTemperatureSection(), "two_temperature", "[temperature] model = " + t);
		if (t != "pump" && t != "probe")
			unused(c.pump_probe == PumpProbeSection(), "pump_probe", "[temperature] model = " + t);
		if (t != "none")
			require(e == "local" || e == "hydrodynamic",
					"[temperature] models combine with [electrons] model = local or hydrodynamic");
		if (t != "none" && e == "hydrodynamic")
			require(t == "two_temperature" && !c.hydrodynamic.split_nonlinear,
					"[electrons] model = hydrodynamic supports [temperature] model = two_temperature only, without split_nonlinear");
		if (c.two_temperature.save_frames && ttm)
			require(t == "two_temperature", "[two_temperature] save_frames needs [temperature] model = two_temperature");
		require(c.two_temperature.rk_stages >= 1, "[two_temperature] rk_stages must be >= 1");
	}
	if (c.output.dft_magnetic)
		require(c.output.dft_xy || c.output.dft_xz || c.output.dft_yz, "[output] dft_magnetic needs at least one of dft_xy, dft_xz, dft_yz");

	auto check_chunks = [&](const std::string& key, const std::string& raw, int nproc, int ncell) {
		if (raw.empty())
			return;
		auto v = parse_int_list("parallel", key, raw);
		require(int(v.size()) == nproc, "[parallel] " + key + " needs exactly " + std::to_string(nproc) + " entries");
		require(!v.empty() && v.front() == 0, "[parallel] " + key + " must start at 0");
		require(std::is_sorted(v.begin(), v.end()) && std::adjacent_find(v.begin(), v.end()) == v.end(),
				"[parallel] " + key + " must be strictly increasing");
		require(v.empty() || v.back() < ncell, "[parallel] " + key + " entries must be < the axis length");
	};
	check_chunks("chunks_x", c.parallel.chunks_x, c.parallel.ranks_x, g.nx + 2 * g.pml_x);
	check_chunks("chunks_y", c.parallel.chunks_y, c.parallel.ranks_y, g.ny + 2 * g.pml_y);
	check_chunks("chunks_z", c.parallel.chunks_z, c.parallel.ranks_z, g.nz + 2 * g.pml_z);
	{
		const char* axes[3] = {"x", "y", "z"};
		const int totals[3] = {g.nx + 2 * g.pml_x, g.ny + 2 * g.pml_y, g.nz + 2 * g.pml_z};
		const int ranks[3] = {c.parallel.ranks_x, c.parallel.ranks_y, c.parallel.ranks_z};
		const std::string* chunk[3] = {&c.parallel.chunks_x, &c.parallel.chunks_y, &c.parallel.chunks_z};
		for (int a = 0; a < 3; ++a)
		{
			try { rank_starts(axes[a], totals[a], ranks[a], *chunk[a], tfsf_planes(c, a)); }
			catch (const validation_error& e) { err.push_back(e.what()); }
		}
	}

	if (process > 1)
	{
		bool any = false;
		for (const auto& d : c.definitions)
			any = any || !(process == 2 ? d.chi2 : d.chi3).empty();
		require(any, std::string("[nonlinear] ") + (process == 2 ? "shg" : "thg") + " needs a material with a " + (process == 2 ? "chi2" : "chi3") +
						 " tensor: chi" + std::to_string(process) + " = <name> in its [material.<name>] section");
	}
	if (lower(c.electrons.model) == "hydrodynamic")
	{
		bool any = false;
		for (const auto& d : c.definitions)
			any = any || d.kind >= materials::Kind::hydrodynamic;
		require(any, "[electrons] model = hydrodynamic needs a material of model hydrodynamic, drude_lorentz or magnetic_lorentz");
	}

	const int sphere = (c.multipoles.sphere_radius > 0) + (c.multipoles.sphere_theta > 0) + (c.multipoles.sphere_phi > 0);
	require(sphere == 0 || sphere == 3, "[multipoles] sphere_radius, sphere_theta and sphere_phi must be set together");

	require(c.probes.size() <= 10, "at most 10 [probe.*] sections are supported");
	for (const auto& p : c.probes)
		require(p.x >= 0 && p.y >= 0 && p.z >= 0, "[probe." + p.name + "] coordinates must be >= 0");

	{
		const std::string shape = lower(c.geometry.shape);
		const ShapeSection* own = find_shape_section(shape);
		for (const auto& sh : SHAPE_SECTIONS)
			if (&sh != own)
				for (const auto& k : sh.keys)
					require(shape_param(c, sh.name, k.key) == 0, std::string("[") + sh.name + "] " + k.key + " is set, but [geometry] shape = " + shape + " does not use that section");
		if (own)
		{
			for (const auto& k : own->keys)
				if (k.pool == &G::shift_y)
					require(c.geometry.shift_y == 0, "[geometry] shift_y is set, but shape = " + shape + " uses [" + shape + "] " + k.key);
		}
	}
	if (lower(c.geometry.shape) == "cylinder_shell")
		require(shape_param(c, "cylinder_shell", "radius_shell") > shape_param(c, "cylinder_shell", "radius_core") && shape_param(c, "cylinder_shell", "radius_core") > 0,
				"[cylinder_shell] needs radius_shell > radius_core > 0");
	if (lower(c.geometry.shape) != "anapole")
		require(c.anapole == AnapoleSection(), "[anapole] has keys set, but [geometry] shape = " + lower(c.geometry.shape) + " does not use that section");
	if (lower(c.geometry.shape) == "anapole")
	{
		const auto& a = c.anapole;
		const double s_gap = (0.5 * a.gap_length / a.patch_edge - a.lobe_offset) / a.lobe_length;
		require(a.patch_edge > 0 && a.thickness > 0, "[anapole] patch_edge and thickness must be > 0");
		require(a.corner_radius >= 0 && 2 * a.corner_radius <= a.patch_edge, "[anapole] corner_radius must be in [0, patch_edge / 2]");
		require(a.lobe_length > 0 && a.lobe_height > 0 && a.lobe_exponent > 0, "[anapole] lobe_length, lobe_height and lobe_exponent must be > 0");
		require(a.gap_length > 0 && a.gap_width > 0, "[anapole] gap_length and gap_width must be > 0");
		require(s_gap > 0 && s_gap < 1, "[anapole] the gap must end on the lobes: 0 < (gap_length / 2 / patch_edge - lobe_offset) / lobe_length < 1");
		require(a.lobe_points >= 8, "[anapole] lobe_points must be >= 8");
		require(a.lobe_offset + a.lobe_length < 0.5, "[anapole] the lobes must end inside the patch: lobe_offset + lobe_length < 1/2");
	}
	for (const auto& m : c.monitors)
		require(m.size_x > 0 && m.size_y > 0 && m.size_z > 0 && m.step > 0,
				"[monitor." + m.name + "] sizes and step must be > 0");
	// source: one type and one signal; the sections of the others must be untouched
	const std::string type = lower(c.source.type);
	const std::string signal = lower(c.source.signal);
	auto only_for = [&](bool untouched, const std::string& section, const std::string& key, const std::string& value) {
		require(untouched, "[" + section + "] has keys set, but [source] " + key + " = " + value + " does not use that section");
	};
	if (type != "plane_wave")
		only_for(c.plane_wave == PlaneWaveSection(), "plane_wave", "type", type);
	if (type != "wave_superposition")
		only_for(c.wave_superposition == WaveSuperpositionSection(), "wave_superposition", "type", type);
	if (type != "gaussian_beam")
		only_for(c.gaussian_beam == GaussianBeamSection(), "gaussian_beam", "type", type);
	if (type != "mode")
		only_for(c.mode == ModeSection(), "mode", "type", type);
	if (type == "mode")
	{
		const auto& m = c.mode;
		require(!m.efield.empty() && !m.hfield.empty(), "[mode] efield and hfield, the MPB field files of the band, are required for [source] type = mode");
		require(m.effective_index > 0, "[mode] effective_index, k over the frequency of the band, is required for [source] type = mode");
		require(m.unit_length > 0, "[mode] unit_length must be > 0");
		require(signal != "raised_cosine" && signal != "raised_cosine_derivative" && signal != "rectangular" && signal != "sinc" && signal != "constant",
				"[source] type = mode needs a signal with a carrier at the mode's wavelength: gaussian, gaussian_alt, continuous_wave or cw_pulse");
		const std::string n = lower(m.propagation_axis);
		const std::string a = lower(m.mpb_x_axis);
		const std::string b = lower(m.mpb_y_axis);
		for (const std::string& t : {a, b})
			require(t.empty() || ((t == "x" || t == "y" || t == "z") && t != n), "[mode] mpb_x_axis and mpb_y_axis must be x, y or z, across the propagation axis");
		require(a.empty() == b.empty() && (a.empty() || a != b), "[mode] set both mpb_x_axis and mpb_y_axis, on different axes, or neither");
	}
	if (type != "plane_wave" && type != "wave_superposition")
		require(c.source.aux_length_factor == 1, "[source] aux_length_factor is for type = plane_wave or wave_superposition");
	auto signal_section = [&](const char* name, bool untouched) {
		if (signal != name)
			only_for(untouched, std::string("signal.") + name, "signal", signal);
	};
	signal_section("raised_cosine", c.signal.raised_cosine == BandLimitedSignal());
	signal_section("raised_cosine_derivative", c.signal.raised_cosine_derivative == BandLimitedSignal());
	signal_section("rectangular", c.signal.rectangular == BandLimitedSignal());
	signal_section("gaussian", c.signal.gaussian == GaussianSignal());
	signal_section("gaussian_alt", c.signal.gaussian_alt == GaussianAltSignal());
	signal_section("continuous_wave", c.signal.continuous_wave == ContinuousWaveSignal());
	signal_section("cw_pulse", c.signal.cw_pulse == CwPulseSignal());
	require(c.source.aux_length_factor >= 1, "[source] aux_length_factor must be >= 1");

	if (type == "plane_wave")
	{
		const auto& p = c.plane_wave;
		const std::string axis = lower(p.polarization_axis);
		require(p.propagation_angle == 0, "[plane_wave] propagation_angle must be 0; oblique incidence is [source] type = wave_superposition");
		require(p.polarization_angle >= -180 && p.polarization_angle <= 180, "[plane_wave] polarization_angle must be in [-180, 180]");
		if (axis == "rhcp" || axis == "lhcp")
		{
			require(p.polarization_angle == 0, "[plane_wave] polarization_angle has no meaning with polarization_axis = " + axis);
			require(signal == "gaussian" || signal == "gaussian_alt" || signal == "continuous_wave" || signal == "cw_pulse",
					"[plane_wave] polarization_axis = " + axis + " needs a signal with a carrier: gaussian, gaussian_alt, continuous_wave or cw_pulse");
		}
	}
	if (type == "dipole")
		require(!c.dipoles.empty(), "[source] type = dipole needs at least one [dipole.<name>] section");
	else
		require(c.dipoles.empty(), "[dipole.*] sections need [source] type = dipole");
	for (const auto& d : c.dipoles)
		require(d.x >= 0 && d.x <= c.grid.nx && d.y >= 0 && d.y <= c.grid.ny && d.z >= 0 && d.z <= c.grid.nz,
				"[dipole." + d.name + "] lies outside the domain");
	for (const auto& b : c.blocks)
	{
		require(b.x_min < b.x_max && b.y_min < b.y_max && b.z_min < b.z_max, "[block." + b.name + "] needs min < max on every axis");
		require(b.count >= 1, "[block." + b.name + "] count must be >= 1");
		require(b.material >= 0, "[block." + b.name + "] material must be >= 0");
	}
	if (type == "wave_superposition")
		require(!c.waves.empty(), "[source] type = wave_superposition needs at least one [wave.<name>] section");
	else
		require(c.waves.empty(), "[wave.*] sections need [source] type = wave_superposition");
	for (const auto& w : c.waves)
	{
		const double n = std::sqrt(w.direction_x * w.direction_x + w.direction_y * w.direction_y + w.direction_z * w.direction_z);
		require(std::fabs(n - 1) < 1e-6, "[wave." + w.name + "] direction must be a unit vector");
	}

	if (type == "gaussian_beam")
	{
		require(c.gaussian_beam.theta_steps > 0 && c.gaussian_beam.phi_steps > 0,
				"[gaussian_beam] theta_steps and phi_steps are required for type = gaussian_beam");
		require(c.gaussian_beam.focal_factor > 0 && c.gaussian_beam.waist_factor > 0 && c.gaussian_beam.theta_max > 0,
				"[gaussian_beam] focal_factor, waist_factor and theta_max must be > 0 for type = gaussian_beam: a zero focal length or waist gives 0/0 amplitudes");
	}

	if (!err.empty())
	{
		std::string msg = "invalid configuration:";
		for (const auto& e : err)
			msg += "\n  " + e;
		throw validation_error(msg);
	}
}

const std::vector<ShapeSection>& shape_sections() { return SHAPE_SECTIONS; }

const ShapeSection* find_shape_section(const std::string& shape)
{
	for (const auto& sh : SHAPE_SECTIONS)
		if (lower(shape) == sh.name)
			return &sh;
	return nullptr;
}

GeometrySection resolved_geometry(const SimConfig& c)
{
	GeometrySection g = c.geometry;
	if (const ShapeSection* own = find_shape_section(c.geometry.shape))
	{
		for (const auto& k : own->keys)
			g.*(k.pool) = shape_param(c, own->name, k.key) + k.offset;
		if (own->finish)
			own->finish(g);
	}
	return g;
}

std::string to_ini(const SimConfig& c)
{
	Writer w;
	visit_sim_config(c, w);
	for (const auto& p : c.probes)
	{
		w.out << "\n[probe." << p.name << "]\nx = " << p.x << "\ny = " << p.y << "\nz = " << p.z << "\n";
	}
	for (auto m : c.monitors)
		visit_monitor(m, w, "monitor." + m.name);
	for (auto p : c.waves)
		visit_wave(p, w, "wave." + p.name);
	for (auto d : c.dipoles)
		visit_dipole(d, w, "dipole." + d.name);
	for (auto b : c.blocks)
		visit_block(b, w, "block." + b.name);
	// the materials, extends resolved, and the tensors they name: the run repeats without the catalog
	std::set<std::string> tensors;
	for (const auto& d : c.definitions)
	{
		if (d.name == "vacuum")
			continue;
		w.out << materials::to_ini(d);
		if (!d.chi2.empty() && tensors.insert("chi2." + d.chi2).second)
			w.out << materials::tensor_to_ini(2, d.chi2, d.chi2_values.data());
		if (!d.chi3.empty() && tensors.insert("chi3." + d.chi3).second)
			w.out << materials::tensor_to_ini(3, d.chi3, d.chi3_values.data());
	}
	return w.out.str();
}

std::string template_ini()
{
	SimConfig c;
	Writer w;
	w.with_docs = true;
	w.out << "; 3D-FDTD configuration template, every key with its default; required keys are marked\n";
	visit_sim_config(c, w);
	w.out << "\n; time-domain probes, one section each, at most 10:\n; [probe.center]\n; x = 100\n; y = 100\n; z = 100\n";
	w.out << "\n; monitor boxes, one section each:\n; [monitor.gap]\n; type = frequency    ; frequency | time\n"
			 "; center_x = 100\n; center_y = 100\n; center_z = 100\n; size_x = 20\n; size_y = 20\n; size_z = 20\n; step = 1\n";
	w.out << "\n; point dipoles of [source] type = dipole, one section each:\n; [dipole.a]\n; x = 100\n; y = 100\n; z = 100\n; e_x = 0\n; e_y = 0\n; e_z = 1\n";
	w.out << "\n; boxes of one material, layers and gratings, drawn before the [geometry] shape in file order:\n; [block.substrate]\n"
			 "; material = glass\n; y_max = -50e-9      ; bounds from the structure centre [m]; omitted: past the domain\n"
			 "; [block.teeth]\n; material = silicon\n; x_min = 0\n; x_max = 500e-9\n; y_min = -400e-9\n; y_max = 0\n; count = 10\n; pitch_x = 1e-6\n";
	w.out << "\n; materials by name: a [material.<name>] section here, or <name>.ini in the catalog, docs/Materials.md.\n"
			 "; The model's keys follow its model; poles are [material.<name>.pole.1], .2, ...\n"
			 "; [material.glass]\n; model = dielectric\n; n = 1.45\n"
			 "; [material.gold_kerr]\n; extends = gold_vial_d2cp\n; chi3 = au_isotropic   ; a [chi3.au_isotropic] section or catalog chi3/au_isotropic.ini\n"
			 "; [chi2.my_crystal]   ; d11 ... d36, d_mn in m/V, chi2 = 2 d\n; d22 = 2.1e-12\n; d31 = 0.15e-12\n";
	w.out << "\n; plane waves of [source] type = wave_superposition, one section each:\n; [wave.a]\n"
			 "; direction_x = 0\n; direction_y = 1\n; direction_z = 0\n; e_x = 0\n; e_y = 0\n; e_z = 1\n; h_x = 1\n; h_y = 0\n; h_z = 0\n";
	return w.out.str();
}

namespace
{
/// the material names of a run in slot order, each once: vacuum, structure, the three media, the
/// shell of the shell shapes, the anapole fills, the blocks
std::vector<std::string> names_in_slot_order(const std::string& shape, std::vector<std::string> named, const std::string& shell,
											 const std::string& aperture, const std::string& gap, const std::vector<std::string>& blocks)
{
	const std::string sh = lower(shape);
	if (sh == "sphere_shell" || sh == "cylinder_shell")
		named.push_back(shell);
	if (sh == "anapole")
		named.push_back(aperture), named.push_back(gap);
	named.insert(named.end(), blocks.begin(), blocks.end());
	std::vector<std::string> out = {"vacuum"};
	for (auto n : named)
	{
		n = lower(n);
		if (!n.empty() && std::find(out.begin(), out.end(), n) == out.end())
			out.push_back(n);
	}
	return out;
}

/// the file of catalog section `sect` in directory `dir`: material.<name> is <name>.ini, chi2.<name> chi2/<name>.ini
std::string catalog_file(const std::string& dir, const std::string& sect)
{
	const auto dot = sect.find('.');
	const std::string family = sect.substr(0, dot);
	const std::string name = sect.substr(dot + 1);
	return family == "material" ? dir + "/" + name + ".ini" : dir + "/" + family + "/" + name + ".ini";
}

std::vector<std::string> split_paths(const std::string& list)
{
	std::vector<std::string> out;
	std::stringstream ss(list);
	std::string item;
	while (std::getline(ss, item, ':'))
		if (!item.empty())
			out.push_back(item);
	return out;
}
} // namespace

std::vector<std::string> material_names(const SimConfig& c)
{
	std::vector<std::string> blocks;
	for (const auto& b : c.blocks)
		blocks.push_back(b.material_name);
	const auto& g = c.geometry;
	return names_in_slot_order(g.shape, {g.material, g.first_medium, g.middle_medium, g.last_medium}, g.shell_material,
							   c.anapole_fill.aperture, c.anapole_fill.gap, blocks);
}

int material_slot(const SimConfig& c, const std::string& name)
{
	const std::string n = lower(name);
	for (size_t s = 0; s < c.definitions.size(); ++s)
		if (c.definitions[s].name == n)
			return int(s);
	return -1;
}

std::vector<std::string> materials_catalog(const io::Ini& ini, const std::string& config_file)
{
	namespace fsys = std::filesystem;
	std::vector<std::string> dirs;
	const fsys::path base = config_file.empty() ? fsys::current_path() : fsys::path(config_file).parent_path();
	for (const auto& d : split_paths(ini.get("materials", "catalog")))
		dirs.push_back(fsys::path(d).is_absolute() ? d : (base / d).string());
	if (const char* env = std::getenv("FDTD_MATERIALS"))
		for (const auto& d : split_paths(env))
			dirs.push_back(d);
	std::error_code ec;
	const fsys::path exe = fsys::read_symlink("/proc/self/exe", ec);
	if (!ec)
		dirs.push_back((exe.parent_path() / "materials").string());
	return dirs;
}

std::string with_catalog(const std::string& ini_text, const std::string& name, const std::vector<std::string>& overrides)
{
	std::string text = ini_text;
	io::Ini ini = build_config(text, name, overrides);
	const std::vector<std::string> dirs = materials_catalog(ini, name);
	std::vector<std::string> blocks;
	for (const auto& sect : ini.sections_in_file_order())
		if (sect.rfind("block.", 0) == 0)
			blocks.push_back(ini.get(sect, "material"));
	auto key = [&](const char* sect, const char* k, const char* dflt) {
		const std::string v = ini.get(sect, k);
		return v.empty() ? std::string(dflt) : v;
	};
	std::vector<std::string> wanted;
	for (const auto& n : names_in_slot_order(key("geometry", "shape", "none"),
											 {key("geometry", "material", "vacuum"), key("geometry", "first_medium", "vacuum"),
											  key("geometry", "middle_medium", "vacuum"), key("geometry", "last_medium", "vacuum")},
											 key("geometry", "shell_material", "vacuum"), key("anapole", "aperture_material", ""),
											 key("anapole", "gap_material", ""), blocks))
		if (n != "vacuum")
			wanted.push_back(materials::material_section(n));
	// every material the file defines is checked, so what it extends and names is read too
	for (const auto& sect : ini.sections())
		if (sect.rfind("material.", 0) == 0 && sect.find('.', 9) == std::string::npos)
			wanted.push_back(sect);
	std::set<std::string> done;
	while (!wanted.empty())
	{
		const std::string sect = wanted.back();
		wanted.pop_back();
		if (!done.insert(sect).second)
			continue;
		io::Ini plain;
		plain.load_string(text, name);
		const bool by_overrides = ini.is_section(sect) && (sect.rfind("material.", 0) != 0 || ini.is_key(sect, "model") || ini.is_key(sect, "extends"));
		if (!plain.is_section(sect))
		{
			std::string found;
			for (const auto& dir : dirs)
				if (std::filesystem::exists(catalog_file(dir, sect)))
				{
					found = catalog_file(dir, sect);
					break;
				}
			if (found.empty() && by_overrides)
				;
			else if (found.empty())
			{
				std::string where;
				for (const auto& dir : dirs)
					where += "\n  " + catalog_file(dir, sect);
				throw std::invalid_argument("[" + sect + "] is neither in the run file nor in the catalog; looked for" + where);
			}
			else
			{
				std::ifstream in(found);
				std::ostringstream ss;
				ss << in.rdbuf();
				text += "\n; " + found + "\n" + ss.str();
				ini = build_config(text, name, overrides);
				io::Ini check;
				check.load_string(text, name);
				if (!check.is_section(sect))
					throw std::invalid_argument(found + " does not define [" + sect + "]");
			}
		}
		if (sect.rfind("material.", 0) == 0)
			for (const auto& r : materials::references(ini, sect.substr(9)))
				wanted.push_back(r);
	}
	return text;
}

} // namespace fdtd::config
