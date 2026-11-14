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
#include "config/exceptions.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <set>
#include <sstream>
#include "materials/codes.h"

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
	{"sphere_shell", {{"radius", &G::width, "[m] radius of the core; the shell ([geometry] shell_material) is 13 cells thicker"}}},
	{"two_spheres", {{"radius", &G::width, "[m]"}, {"separation", &G::shift_y, "[m] centre to centre, along y"}}},
	{"sphere_dimer", {{"radius", &G::width, "[m]"}, {"spacing", &G::length, "[m] centre to centre, along z"}}},
	{"sphere_trimer", {{"radius", &G::width, "[m]"}, {"spacing", &G::length, "[m] centre to centre, along z"}}},
	{"sphere_quadrumer", {{"radius", &G::width, "[m]"}, {"spacing", &G::length, "[m] centre to centre, along z"}}},
	{"cylinder", {{"radius", &G::width, "[m] cylinder along z"}, {"length", &G::length, "[m] along z; 0 = through the whole domain"}}},
	{"cylinder_3d", {{"radius", &G::width, "[m] disc in the x-z plane"}, {"height", &G::length, "[m] along y"}}},
	{"cylinder_shell", {{"radius_core", &G::width, "[m] radius of the core; both cylinders are along z, on the same axis"},
						{"core_length", &G::length, "[m] of the core along z; 0 = through the whole domain"},
						{"radius_shell", &G::length_2, "[m] outer radius of the shell ([geometry] shell_material); the core is drawn over it"},
						{"shell_length", &G::arm_length, "[m] of the shell along z; 0 = through the whole domain; shorter than the core leaves the core ends bare"}}},
	{"ellipse", {{"radius_x", &G::width, "[m] semi-axis along x"}, {"radius_z", &G::length, "[m] semi-axis along z"}, {"thickness", &G::thickness, "[m] along y"},
				 {"rotation", &G::rotation, "[deg] positive turns counter-clockwise in the top view (x right, z up)"}}},
	{"hole", {{"radius", &G::width, "[m] of the hole"}, {"thickness", &G::thickness, "[m] of the film that fills the x-z plane, along y"}}},
	{"dipole", {{"width", &G::width, "[m] arm width, along x"}, {"thickness", &G::thickness, "[m] along y"}, {"gap", &G::gap, "[m] between the arms, along z"},
				{"top_end", &G::length, "[m] end of the top arm, measured from the centre of the gap (+z)"},
				{"bottom_end", &G::length_2, "[m] end of the bottom arm, measured from the centre of the gap (-z)"}}},
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
					  {"corner_radius", &G::corner_radius, "[m] of every edge"}, {"rotation", &G::rotation, "[deg] 0: arms along z; positive turns counter-clockwise in the top view (x right, z up)", ROUND_DIPOLE_ZERO}}},
	{"one_side_round_dipole", {{"width", &G::width, "[m] arm width; the outer ends are half discs of this diameter"}, {"thickness", &G::thickness, "[m] along y"},
							   {"gap", &G::gap, "[m] between the arms"},
							   {"top_end", &G::length, "[m] centre of the top end disc, measured from the centre of the gap"},
							   {"bottom_end", &G::arm_length, "[m] centre of the bottom end disc, measured from the centre of the gap"},
							   {"rotation", &G::rotation, "[deg] 0: arms along z; positive turns counter-clockwise in the top view (x right, z up)", ROUND_DIPOLE_ZERO}}},
	{"round_hybrid", {{"width", &G::width, "[m] arm width"}, {"thickness", &G::thickness, "[m] along y"}, {"gap", &G::gap, "[m] between the arms (filled with the gap material)"},
					  {"arm_length", &G::length_2, "[m] length of each arm"}, {"corner_radius", &G::corner_radius, "[m] of every edge"},
					  {"rotation", &G::rotation, "[deg] 0: arms along z; positive turns counter-clockwise in the top view (x right, z up)", ROUND_DIPOLE_ZERO}},
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
				   {"length_x", &G::length, "[m] of the bar along x (at the -z end), from the outer corner to the centre of its round end"},
				   {"length_z", &G::length_2, "[m] of the bar along z (at the -x end), from the outer corner to the centre of its round end"}}},
	{"lateral_dipole", {{"width", &G::width, "[m] width of each bar, along x"}, {"thickness", &G::thickness, "[m] along y"},
						{"gap", &G::gap, "[m] between the two bars, along x"},
						{"left_length", &G::length, "[m] of the bar at -x; it reaches toward +z"},
						{"right_length", &G::length_2, "[m] of the bar at +x; it reaches toward -z"},
						{"overlap", &G::arm_length, "[m] length along z over which the two bars face each other; the facing corners are cut at 45 deg"}}},
	{"round_butterfly", {{"width", &G::width, "[m] bar width"}, {"thickness", &G::thickness, "[m] along y"}, {"gap", &G::gap, "[m] between the two halves; each half is one z bar and one x bar joined by a diagonal neck"},
						 {"length_z", &G::length_2, "[m] of each bar along z, measured from width / 2 off the centre"}, {"length_x", &G::arm_width, "[m] of each bar along x, measured the same way"},
						 {"corner_radius", &G::corner_radius, "[m] of every edge and of the two arcs"}, {"rotation", &G::rotation, "[deg] 0: bars along x and z; positive turns counter-clockwise in the top view (x right, z up)", BUTTERFLY_ZERO}}},
	{"round_tripole", {{"width", &G::width, "[m] bar width"}, {"thickness", &G::thickness, "[m] along y"},
					   {"center_length", &G::length, "[m] of the middle bar; the three bars are in line"}, {"top_length", &G::length_2, "[m] of the top bar (+z at rotation 0)"}, {"bottom_length", &G::arm_width, "[m] of the bottom bar"},
					   {"top_gap", &G::arm_length, "[m] between the middle and the top bar"}, {"bottom_gap", &G::gap, "[m] between the middle and the bottom bar"},
					   {"corner_radius", &G::corner_radius, "[m] of every edge"}, {"rotation", &G::rotation, "[deg] 0: bars along z; positive turns counter-clockwise in the top view (x right, z up)", ROUND_DIPOLE_ZERO}}},
	{"x_shape", {{"thickness", &G::thickness, "[m] along y"},
				 {"fixed_length", &G::length_2, "[m] of the fixed bar, along x"}, {"fixed_width", &G::width, "[m] of the fixed bar"},
				 {"rotated_length", &G::length, "[m] of the second bar, which crosses the fixed one at its centre"}, {"rotated_width", &G::gap, "[m] of the second bar"},
				 {"corner_radius", &G::corner_radius, "[m] of every edge"},
				 {"angle", &G::rotation, "[deg] from the fixed bar to the second bar, counter-clockwise in the top view (x right, z up); 90 is a cross"}}},
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
const std::vector<std::string> HARMONIC_MODES = {"none", "antenna", "bulk"};

/** Every scalar setting, exactly once. V must provide
 *    field(section, key, T& ref, doc, required)
 *    choice(section, key, std::string& ref, choices, doc)
 *    refractive_indices(map) */
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
	v.field(s, "e_x", d.e_x, "amplitude of the injected Ex (times the pulse)");
	v.field(s, "e_y", d.e_y, "");
	v.field(s, "e_z", d.e_z, "");
}

/// [block.<name>]: one box of a layered or periodic structure
template <class B, class V>
void visit_block(B& b, V& v, const std::string& sect)
{
	const char* s = sect.c_str();
	v.field(s, "material", b.material, "material code", true);
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
	v.field(s, "h_x", w.h_x, "H amplitude (vacuum units; scaled by n of the first medium)");
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

	v.field("grid", "nx", c.grid.nx, "interior cells along x (PML added on top)", true);
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
	v.field("grid", "staircase", c.grid.staircase, "staircased (true) or per-component (false) material assignment");

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
	v.field("source", "allow_truncated_signal", c.source.allow_truncated_signal, "run even when the time steps do not cover the signal (the run stops otherwise)");
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
	v.field("mode", "efield", c.mode.efield, "MPB E-field file of the band (e.k01.b01.h5), relative to [paths] input_dir");
	v.field("mode", "hfield", c.mode.hfield, "MPB H-field file of the band (h.k01.b01.h5)");
	v.field("mode", "unit_length", c.mode.unit_length, "[m] MPB's length unit a");
	v.field("mode", "effective_index", c.mode.effective_index, "n_eff = k / frequency of the band (MPB units)");
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
	v.choice("spectrum", "spacing", c.spectrum.spacing, {"wavelength", "frequency"}, "uniform step in");
	v.choice("spectrum", "harmonic_band", c.spectrum.harmonic_band, {"none", "second", "third"}, "extra DFT band at 2w / 3w of every point");

	v.field("output", "dft_xy", c.output.dft_xy, "frequency-domain E on the XY plot plane (fields/E_<wl>nm_XY.h5)");
	v.field("output", "dft_xz", c.output.dft_xz, "frequency-domain E on the XZ plot plane");
	v.field("output", "dft_yz", c.output.dft_yz, "frequency-domain E on the YZ plot plane");
	v.field("output", "dft_magnetic", c.output.dft_magnetic, "also H on the enabled planes (fields/H_<wl>nm_<plane>.h5)");
	v.field("output", "volume_dft", c.output.volume_dft, "full-volume E/H DFT");
	v.field("output", "frames_xy", c.output.frames_xy, "time-domain XY frame period [steps]; 0 = off");
	v.field("output", "frames_xz", c.output.frames_xz, "");
	v.field("output", "frames_yz", c.output.frames_yz, "");
	v.field("output", "frames_3d", c.output.frames_3d, "");
	v.field("output", "phase_normalize", c.output.phase_normalize, "normalize the DFT phase to the source");
	v.field("output", "xdmf", c.output.xdmf, "write a .xmf next to every HDF5 output, for ParaView / VisIt");
	v.field("output", "sigma_absorption", c.output.sigma_absorption, "absorption from the sigma|E|^2 volume integral");

	v.field("multipoles", "exact", c.multipoles.exact, "exact multipole decomposition (Alaee et al. 2018)");
	v.field("multipoles", "lwa", c.multipoles.lwa, "long-wavelength multipoles (kr -> 0 limit of the same integrals); any material model");
	v.field("multipoles", "sphere_radius", c.multipoles.sphere_radius, "legacy sphere-sampled multipoles: radius [cells]; 0 = off");
	v.field("multipoles", "sphere_theta", c.multipoles.sphere_theta, "");
	v.field("multipoles", "sphere_phi", c.multipoles.sphere_phi, "");

	v.choice("geometry", "shape", c.geometry.shape, shape_names(), "structure");
	v.field("geometry", "material", c.geometry.material, "material of the structure");
	v.field("geometry", "first_medium", c.geometry.first_medium, "entry layer material");
	v.field("geometry", "middle_medium", c.geometry.middle_medium, "middle layer material");
	v.field("geometry", "last_medium", c.geometry.last_medium, "exit layer material");
	v.field("geometry", "middle_thickness", c.geometry.middle_thickness, "[m]");
	v.field("geometry", "shift_y", c.geometry.shift_y, "structure offset along y [m]");
	v.field("geometry", "shell_material", c.geometry.shell_material, "");

	for (const auto& sh : SHAPE_SECTIONS)
		for (const auto& k : sh.keys)
			v.field(sh.name, k.key, shape_param(c, sh.name, k.key), k.doc);

	v.field("materials", "carrier_density", c.materials.carrier_density, "");
	v.field("materials", "carrier_density_back", c.materials.carrier_density_back, "");
	v.field("materials", "stretch_x", c.materials.stretch_x, "graphene (material 102): scales the plasma frequency and critical-point amplitudes");
	v.refractive_indices(c.materials.refractive_index);

	v.field("custom_material", "eps_inf", c.custom_material.eps_inf, "");
	v.field("custom_material", "plasma_frequency", c.custom_material.plasma_frequency, "");
	v.field("custom_material", "plasma_damping", c.custom_material.plasma_damping, "");
	v.field("custom_material", "cp_frequency", c.custom_material.cp_frequency, "critical point 1");
	v.field("custom_material", "cp_frequency_2", c.custom_material.cp_frequency_2, "critical point 2");
	v.field("custom_material", "cp_damping", c.custom_material.cp_damping, "");
	v.field("custom_material", "cp_damping_2", c.custom_material.cp_damping_2, "");
	v.field("custom_material", "cp_amplitude", c.custom_material.cp_amplitude, "");
	v.field("custom_material", "cp_amplitude_2", c.custom_material.cp_amplitude_2, "");
	v.field("custom_material", "cp_phase", c.custom_material.cp_phase, "");
	v.field("custom_material", "cp_phase_2", c.custom_material.cp_phase_2, "");
	v.field("custom_material", "fermi_velocity", c.custom_material.fermi_velocity, "");
	v.field("custom_material", "electron_density", c.custom_material.electron_density, "");

	v.field("anapole", "patch_edge", c.anapole.patch_edge, "[m] edge of the square patch (shape = anapole)");
	v.field("anapole", "thickness", c.anapole.thickness, "[m] film thickness along y");
	v.field("anapole", "corner_radius", c.anapole.corner_radius, "[m] fillet radius of the patch corners");
	v.field("anapole", "lobe_a", c.anapole.lobe_a, "lobe length / patch_edge");
	v.field("anapole", "lobe_b", c.anapole.lobe_b, "lobe height / patch_edge");
	v.field("anapole", "lobe_p", c.anapole.lobe_p, "lobe exponent");
	v.field("anapole", "lobe_x0", c.anapole.lobe_x0, "lobe start / patch_edge");
	v.field("anapole", "neck_half_length", c.anapole.neck_half_length, "[m] neck half-length along the lobe axis");
	v.field("anapole", "neck_half_height", c.anapole.neck_half_height, "[m] neck half-height across the lobe axis");
	v.field("anapole", "rotation", c.anapole.rotation, "[deg] in-plane rotation (x-z)");
	v.field("anapole", "lobe_material", c.anapole.lobe_material, "fill of the lobes; -1 = background");
	v.field("anapole", "neck_material", c.anapole.neck_material, "fill of the neck; -1 = background");
	v.field("anapole", "dimer_material", c.anapole.dimer_material, "two cylinders in the neck; -1 = none");
	v.field("anapole", "dimer_radius", c.anapole.dimer_radius, "[m]");
	v.field("anapole", "dimer_gap", c.anapole.dimer_gap, "[m]");
	v.field("anapole", "lobe_points", c.anapole.lobe_points, "polygon vertices per lobe edge");

	v.choice("nonlinear", "shg", c.nonlinear.shg, HARMONIC_MODES, "second-harmonic generation: in the structure (antenna) or the layers (bulk)");
	v.choice("nonlinear", "thg", c.nonlinear.thg, HARMONIC_MODES, "third-harmonic generation: in the structure (antenna) or the layers (bulk)");
	v.choice("nonlinear", "update", c.nonlinear.update, {"field", "polarization"}, "how the nonlinearity enters the update");
	v.field("nonlinear", "chi_structure", c.nonlinear.chi_structure, "tensor of the structure: a [chi2|chi3.<name>] section or a preset; unset = zero");
	v.field("nonlinear", "chi_first", c.nonlinear.chi_first, "tensor name of the first medium");
	v.field("nonlinear", "chi_middle", c.nonlinear.chi_middle, "tensor name of the middle medium");
	v.field("nonlinear", "chi_last", c.nonlinear.chi_last, "tensor name of the last medium");
	v.field("nonlinear", "graphene_optimization", c.nonlinear.graphene_optimization, "");

	v.choice("electrons", "model", c.electrons.model, {"local", "nonlocal", "gnor", "hydrodynamic"}, "free-electron response; the hydrodynamic keys are in [hydrodynamic]");
	v.field("electrons", "diffusion", c.electrons.diffusion, "diffusion coefficient (nonlocal, gnor, hydrodynamic)");
	v.field("electrons", "fermi_velocity", c.electrons.fermi_velocity, "(nonlocal, gnor, hydrodynamic)");

	v.field("hydrodynamic", "pressure", c.hydrodynamic.pressure, "quantum pressure term");
	v.field("hydrodynamic", "convection", c.hydrodynamic.convection, "convective term");
	v.field("hydrodynamic", "magnetic", c.hydrodynamic.magnetic, "magnetic Lorentz-force term");
	v.field("hydrodynamic", "split_nonlinear", c.hydrodynamic.split_nonlinear, "separate linear / nonlinear DFTs");

	v.choice("temperature", "model", c.temperature.model, {"none", "two_temperature", "load", "pump", "probe"},
			 "electron temperature: two-temperature model, load <input_dir>/electron_temperature.h5, pump-probe");

	v.field("two_temperature", "save_frames", c.two_temperature.save_frames, "write temperature frames (model = two_temperature)");
	v.field("two_temperature", "rk_stages", c.two_temperature.rk_stages, "Runge-Kutta stages");

	v.field("pump_probe", "delay", c.pump_probe.delay, "[s]");
}

// Value codecs

void parse_value(const Config& ini, const std::string& s, const std::string& k, int& out)
{
	out = Config::to_integer(s, k, ini.get(s, k));
}
void parse_value(const Config& ini, const std::string& s, const std::string& k, double& out)
{
	out = Config::to_double(s, k, ini.get(s, k));
}
void parse_value(const Config& ini, const std::string& s, const std::string& k, bool& out)
{
	out = ini.get_bool(s, k);
}
void parse_value(const Config& ini, const std::string& s, const std::string& k, std::string& out)
{
	out = Config::unquote(ini.get(s, k));
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

std::string lower(std::string s)
{
	for (auto& ch : s)
		ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
	return s;
}

// Visitors

struct Loader
{
	const Config& ini;
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

	void refractive_indices(std::map<int, double>&) {} // loaded separately, dynamic keys

	static void check_choice(const char* s, const char* k, const std::string& ref, const std::vector<std::string>& choices)
	{
		if (std::find(choices.begin(), choices.end(), ref) == choices.end())
		{
			std::string all;
			for (const auto& c : choices)
				all += (all.empty() ? "" : ", ") + c;
			throw syntax_error(s, k, "invalid value '" + ref + "' (expected one of: " + all + ")");
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
			out << "; " << (required ? "(required) " : "") << (doc ? doc : "") << "\n";
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

	void refractive_indices(const std::map<int, double>& n)
	{
		header("materials");
		if (with_docs)
			out << "; refractive index of dispersionless material i (i = 0..99): n<i> = value, e.g.\n; n4 = 1.45\n";
		for (const auto& [idx, val] : n)
			out << "n" << idx << " = " << format_value(val) << "\n";
	}
};

/// "n12" -> 12 ; anything else -> -1
int refractive_index_key(const std::string& key)
{
	if (key.size() < 2 || key[0] != 'n')
		return -1;
	for (size_t i = 1; i < key.size(); ++i)
		if (!std::isdigit(static_cast<unsigned char>(key[i])))
			return -1;
	return std::atoi(key.c_str() + 1);
}

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
									" cells along " + axis + " (PML included); use at most " + std::to_string(total / 2));
	std::vector<int> start;
	if (!chunks.empty())
	{
		start = parse_int_list("parallel", key, chunks);
		for (size_t r = 1; r < start.size(); ++r)
			if (on_plane(start[r]))
				throw validation_error("[parallel] " + key + ": the rank boundary at cell " + std::to_string(start[r]) +
											" lies on the TF/SF plane along " + axis +
											"; move it by one cell (e.g. " + std::to_string(start[r] + 1) + ")");
		return start;
	}
	const int base = total / ranks, extra = total % ranks;
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

int harmonic_mode(const std::string& name)
{
	return lookup({{"none", 0}, {"antenna", 1}, {"bulk", 2}}, name, "[nonlinear] shg / thg");
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
	if (lower(n.shg) != "none") return 2;
	if (lower(n.thg) != "none") return 3;
	return 1;
}

std::vector<double> spectrum_omegas(const SpectrumSection& s, double c0, double two_pi)
{
	std::vector<double> w;
	const int n = s.points;
	for (int m = 0; m < n; ++m)
	{
		if (n == 1)
			w.push_back(two_pi * (c0 / s.lambda_min));
		else if (lower(s.spacing) == "frequency")
		{
			const double f0 = c0 / s.lambda_max, f1 = c0 / s.lambda_min;
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
		out.push_back(Config::to_integer(sect, key, tok));
	return out;
}

bool split_override(const std::string& spec, std::string& section, std::string& key, std::string& value)
{
	const auto eq = spec.find('=');
	if (eq == std::string::npos)
		return false;
	std::string path = spec.substr(0, eq);
	value = spec.substr(eq + 1);
	Config::trim(path);
	Config::trim(value);
	// split at the LAST dot so that dynamic sections, "probe.top.x", work
	const auto dot = path.rfind('.');
	if (dot == std::string::npos || dot == 0 || dot + 1 == path.size())
		return false;
	section = path.substr(0, dot);
	key = path.substr(dot + 1);
	return true;
}

Config build_config(const std::string& ini_text, const std::string& name, const std::vector<std::string>& overrides,
						const std::string& input_dir, const std::string& output_dir)
{
	Config ini;
	ini.load_string(ini_text, name.empty() ? "<config>" : name);
	for (const auto& spec : overrides)
	{
		std::string s, k, v;
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

SimConfig load_sim_config(const Config& ini)
{
	SimConfig c;
	Loader loader{ini, {}, {}};
	visit_sim_config(c, loader);

	// dynamic keys
	for (const auto& key : ini.keys("materials"))
	{
		const int idx = refractive_index_key(key);
		if (idx < 0)
			continue;
		if (idx > 99)
			throw syntax_error("materials", key, "refractive indices are defined for material 0..99");
		c.materials.refractive_index[idx] = Config::to_double("materials", key, ini.get("materials", key));
	}

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
		else if (sect.rfind("chi2.", 0) == 0 || sect.rfind("chi3.", 0) == 0)
		{
			ChiTensor t;
			t.order = sect[3] - '0';
			t.name = sect.substr(5);
			for (const auto& key : ini.keys(sect))
			{
				const int i = chi_component_index(t.order, key);
				if (i < 0)
					throw syntax_error(sect, key, t.order == 2 ? "chi2 components are d11 ... d36 (d_mn in m/V)"
																   : "chi3 components are xxxx, xxyy, ... zyzy (in m^2/V^2)");
				t.values[i] = Config::to_double(sect, key, ini.get(sect, key));
			}
			c.chi_tensors.push_back(t);
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
		if (sect.rfind("chi2.", 0) == 0 || sect.rfind("chi3.", 0) == 0)
			continue; // components checked above
		for (const auto& key : ini.keys(sect))
		{
			if (is_probe && (key == "x" || key == "y" || key == "z"))
				continue;
			if (sect == "materials" && refractive_index_key(key) >= 0)
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
			"[boxes] scattering_* must be >= 0 (tfsf_* too small?)");

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
	require(lower(c.nonlinear.shg) == "none" || lower(c.nonlinear.thg) == "none", "[nonlinear] shg and thg cannot both be on");
	if (process <= 1)
		require(c.nonlinear.chi_structure.empty() && c.nonlinear.chi_first.empty() && c.nonlinear.chi_middle.empty() && c.nonlinear.chi_last.empty(),
				"[nonlinear] chi_* are set, but shg and thg are none");
	if (band > 1 && lower(c.electrons.model) == "local" && lower(c.temperature.model) == "none")
		require(process == band, "[spectrum] harmonic_band = " + s.harmonic_band + " requires [nonlinear] " + (band == 2 ? "shg" : "thg") + " = antenna or bulk");
	if (process > 1)
		require(band == process, std::string("[nonlinear] ") + (process == 2 ? "shg" : "thg") +
									 " computes the harmonic: set [spectrum] harmonic_band = " + (process == 2 ? "second" : "third"));
	{
		const std::string e = lower(c.electrons.model), t = lower(c.temperature.model);
		auto unused = [&](bool untouched, const char* section, const std::string& selector) {
			require(untouched, std::string("[") + section + "] has keys set, but " + selector + " does not use that section");
		};
		if (e == "local")
			require(c.electrons.diffusion == ElectronsSection().diffusion && c.electrons.fermi_velocity == ElectronsSection().fermi_velocity,
					"[electrons] diffusion and fermi_velocity need model = nonlocal, gnor or hydrodynamic");
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
	// the nonlinear energy balance, Parseval, needs a constant frequency step in a harmonic band
	if (harmonic_order(s.harmonic_band) > 1 && s.points > 1)
		require(lower(s.spacing) == "frequency", "[spectrum] a harmonic band over several points needs spacing = frequency");

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

	for (int mat : {c.geometry.material, c.geometry.first_medium, c.geometry.middle_medium, c.geometry.last_medium})
		if (mat > 0 && mat < 100)
			require(c.materials.refractive_index.count(mat) > 0,
					"[materials] n" + std::to_string(mat) + " is required (material " + std::to_string(mat) + " is in use)");

	// [custom_material] fills one material code: the user hydrodynamic metal, slot 2 after the
	// modified-Lorentz range, src/materials/parameters.cc
	{
		constexpr int custom = fdtd::materials::DispMateOffsetD + fdtd::materials::ModLorentz + 2;
		bool used = false;
		for (int mat : {c.geometry.material, c.geometry.first_medium, c.geometry.middle_medium, c.geometry.last_medium, c.geometry.shell_material,
						c.anapole.lobe_material, c.anapole.neck_material, c.anapole.dimer_material})
			used = used || mat == custom;
		for (const auto &b : c.blocks)
			used = used || b.material == custom;
		require(used || c.custom_material == CustomMaterialSection(),
				"[custom_material] has keys set, but no structure, layer, block or anapole part uses its material code " + std::to_string(custom));
	}

	if (process > 1)
	{
		const auto& n = c.nonlinear;
		require(!n.chi_structure.empty() || !n.chi_first.empty() || !n.chi_middle.empty() || !n.chi_last.empty(),
				"[nonlinear] SHG/THG needs a tensor: set chi_structure (and chi_first, chi_middle, chi_last for the media; unset = zero)");
		const char* family = process == 2 ? "chi2" : "chi3";
		for (const auto& [key, name] : {std::pair<const char*, const std::string&>{"chi_structure", n.chi_structure},
										{"chi_first", n.chi_first}, {"chi_middle", n.chi_middle}, {"chi_last", n.chi_last}})
		{
			if (name.empty() || find_chi_tensor(process, name, c.chi_tensors))
				continue;
			std::string presets;
			for (const auto& p : chi_preset_names(process))
				presets += (presets.empty() ? "" : ", ") + p;
			require(false, std::string("[nonlinear] ") + key + " = " + name + ": no [" + family + "." + name +
							   "] section and no preset of that name (" + family + " presets: " + presets + ")");
		}
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
		const double s_neck = (a.neck_half_length / a.patch_edge - a.lobe_x0) / a.lobe_a;
		require(a.patch_edge > 0 && a.thickness > 0, "[anapole] patch_edge and thickness must be > 0");
		require(a.corner_radius >= 0 && 2 * a.corner_radius <= a.patch_edge, "[anapole] corner_radius must be in [0, patch_edge / 2]");
		require(a.lobe_a > 0 && a.lobe_b > 0 && a.lobe_p > 0, "[anapole] lobe_a, lobe_b and lobe_p must be > 0");
		require(s_neck > 0 && s_neck < 1, "[anapole] the neck face |u| = neck_half_length must cut the lobes (0 < (neck_half_length / patch_edge - lobe_x0) / lobe_a < 1)");
		require(a.neck_half_height > 0 && a.lobe_points >= 8, "[anapole] neck_half_height must be > 0 and lobe_points >= 8");
		require((a.lobe_x0 + a.lobe_a) * a.patch_edge < 0.5 * a.patch_edge, "[anapole] the lobes must end inside the patch (lobe_x0 + lobe_a < 1/2)");
	}
	for (const auto& m : c.monitors)
		require(m.size_x > 0 && m.size_y > 0 && m.size_z > 0 && m.step > 0,
				"[monitor." + m.name + "] sizes and step must be > 0");
	// source: one type and one signal; the sections of the others must be untouched
	const std::string type = lower(c.source.type), signal = lower(c.source.signal);
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
		require(!m.efield.empty() && !m.hfield.empty(), "[mode] efield and hfield (MPB field files of the band) are required for [source] type = mode");
		require(m.effective_index > 0, "[mode] effective_index (k / frequency of the band) is required for [source] type = mode");
		require(m.unit_length > 0, "[mode] unit_length must be > 0");
		require(signal != "raised_cosine" && signal != "raised_cosine_derivative" && signal != "rectangular" && signal != "sinc" && signal != "constant",
				"[source] type = mode needs a signal with a carrier at the mode's wavelength (gaussian, gaussian_alt, continuous_wave, cw_pulse)");
		const std::string n = lower(m.propagation_axis), a = lower(m.mpb_x_axis), b = lower(m.mpb_y_axis);
		for (const std::string& t : {a, b})
			require(t.empty() || ((t == "x" || t == "y" || t == "z") && t != n), "[mode] mpb_x_axis and mpb_y_axis must be x, y or z, across the propagation axis");
		require(a.empty() == b.empty() && (a.empty() || a != b), "[mode] set both mpb_x_axis and mpb_y_axis (different axes) or neither");
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
		require(p.propagation_angle == 0, "[plane_wave] propagation_angle must be 0 (oblique incidence: [source] type = wave_superposition)");
		require(p.polarization_angle >= -180 && p.polarization_angle <= 180, "[plane_wave] polarization_angle must be in [-180, 180]");
		if (axis == "rhcp" || axis == "lhcp")
		{
			require(p.polarization_angle == 0, "[plane_wave] polarization_angle has no meaning with polarization_axis = " + axis);
			require(signal == "gaussian" || signal == "gaussian_alt" || signal == "continuous_wave" || signal == "cw_pulse",
					"[plane_wave] polarization_axis = " + axis + " needs a signal with a carrier (gaussian, gaussian_alt, continuous_wave, cw_pulse)");
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
				"[gaussian_beam] focal_factor, waist_factor and theta_max must be > 0 for type = gaussian_beam (a zero focal length or waist gives 0/0 amplitudes)");
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
	for (const auto& t : c.chi_tensors)
	{
		w.out << "\n[chi" << t.order << "." << t.name << "]\n";
		const auto& names = chi_components(t.order);
		for (size_t i = 0; i < names.size(); ++i)
			if (t.values[i] != 0)
				w.out << names[i] << " = " << format_value(t.values[i]) << "\n";
	}
	return w.out.str();
}

static std::string join_names(const std::vector<std::string>& v)
{
	std::string s;
	for (const auto& n : v)
		s += (s.empty() ? "" : ", ") + n;
	return s;
}

std::string template_ini()
{
	SimConfig c;
	Writer w;
	w.with_docs = true;
	w.out << "; 3D-FDTD configuration template (generated by `fdtd template`).\n"
			 "; Every key is optional unless marked (required); values shown are the defaults.\n"
			 "; Any key can be overridden from the command line: --set section.key=value\n";
	visit_sim_config(c, w);
	w.out << "\n; time-domain probes, one section each (at most 10):\n; [probe.center]\n; x = 100\n; y = 100\n; z = 100\n";
	w.out << "\n; monitor boxes, one section each:\n; [monitor.gap]\n; type = frequency    ; frequency | time\n"
			 "; center_x = 100\n; center_y = 100\n; center_z = 100\n; size_x = 20\n; size_y = 20\n; size_z = 20\n; step = 1\n";
	w.out << "\n; point dipoles of [source] type = dipole, one section each:\n; [dipole.a]\n; x = 100\n; y = 100\n; z = 100\n; e_x = 0\n; e_y = 0\n; e_z = 1\n";
	w.out << "\n; boxes of one material, drawn before the [geometry] shape in file order (layers, gratings):\n; [block.substrate]\n"
			 "; material = 1\n; y_max = -50e-9      ; bounds from the structure centre [m]; omitted: past the domain\n"
			 "; [block.teeth]\n; material = 2\n; x_min = 0\n; x_max = 500e-9\n; y_min = -400e-9\n; y_max = 0\n; count = 10\n; pitch_x = 1e-6\n";
	w.out << "\n; nonlinear tensors named by [nonlinear] chi_structure / chi_first / chi_middle / chi_last, one section each;\n"
			 "; components not listed are zero. Built-in presets need no section (a section of the same name replaces one):\n"
			 "; chi2: " + join_names(chi_preset_names(2)) + "; chi3: " + join_names(chi_preset_names(3)) + "\n"
			 "; [chi2.my_crystal]   ; d11 ... d36 (d_mn in m/V; chi2 = 2 d)\n; d22 = 2.1e-12\n; d31 = 0.15e-12\n"
			 "; [chi3.my_metal]     ; xxxx xxyy xyyx xyxy xxzz xzzx xzxz yyyy yyxx yxxy yxyx yyzz yzzy yzyz zzzz zzxx zxxz zxzx zzyy zyyz zyzy (m^2/V^2)\n"
			 "; xxxx = 7.6e-19\n; yyyy = 7.6e-19\n; zzzz = 7.6e-19\n";
	w.out << "\n; plane waves of [source] type = wave_superposition, one section each:\n; [wave.a]\n"
			 "; direction_x = 0\n; direction_y = 1\n; direction_z = 0\n; e_x = 0\n; e_y = 0\n; e_z = 1\n; h_x = 1\n; h_y = 0\n; h_z = 0\n";
	return w.out.str();
}

} // namespace fdtd::config
