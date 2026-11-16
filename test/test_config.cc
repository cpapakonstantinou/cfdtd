#include <algorithm>
// Unit tests for fdtd::config::Config, fdtd::config::SimConfig and the CLI layer, no MPI needed.
#include "config/exceptions.h"
#include "config/sim_config.h"

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>

static int failures = 0;
#define CHECK(cond)                                                              \
	do {                                                                         \
		if (!(cond)) { std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); ++failures; } \
	} while (0)

template <class E>
static bool throws(const std::function<void()>& f, const std::string& needle = "")
{
	try { f(); }
	catch (const E& e)
	{
		if (needle.empty() || std::string(e.what()).find(needle) != std::string::npos)
			return true;
		std::printf("  wrong message: %s\n", e.what());
		return false;
	}
	catch (const std::exception& e) { std::printf("  wrong exception: %s\n", e.what()); return false; }
	return false;
}

static const char* BASIC = R"INI(
; minimal run
[grid]
nx = 200
NY = 100          ; keys are case-insensitive
nz = 50
dx = 2e-9
dy = 2e-9
dz = 4e-9
periodic_x = yes

[parallel]
ranks_x = 2

[spectrum]
lambda_min = 500e-9
lambda_max = 900e-9
points = 41

[source]
signal = Gaussian

[signal.gaussian]
center_wavelength = 800e-9

[materials]
n4 = 1.45

[geometry]
material = 4

[multipoles]
exact = true

[probe.center]
x = 10
y = 20
z = 30
)INI";

/// load BASIC-like text with overrides given as {"--set", "section.key=value", ...}
static fdtd::config::SimConfig load(const std::string& text, std::vector<const char*> args = {})
{
	std::vector<std::string> overrides;
	for (std::size_t n = 1; n < args.size(); n += 2)
		overrides.push_back(args[n]);
	auto c = fdtd::config::load_sim_config(fdtd::config::build_config(text, "x.ini", overrides));
	fdtd::config::validate_sim_config(c);
	return c;
}

int main()
{
	// basic load + derived defaults
	auto c = load(BASIC);
	CHECK(c.grid.nx == 200 && c.grid.ny == 100 && c.grid.nz == 50);
	CHECK(c.grid.dx == 2e-9 && c.grid.dz == 4e-9);
	CHECK(c.parallel.ranks_x == 2 && c.parallel.ranks_y == 1);
	CHECK(c.grid.periodic_x && !c.grid.periodic_z);
	CHECK(c.grid.origin_x == 100 && c.grid.origin_y == 50 && c.grid.origin_z == 25);
	CHECK(c.boxes.scattering_x == 6 && c.boxes.absorption_z == 16);
	CHECK(c.source.signal == "gaussian" && fdtd::config::pulse_code(c.source.signal) == 6);
	CHECK(c.signal.gaussian.center_wavelength == 800e-9 && c.plane_wave.polarization_axis == "z");
	CHECK(c.source.aux_steps == c.grid.time_steps);
	CHECK(c.materials.refractive_index.at(4) == 1.45);
	CHECK(c.multipoles.exact);
	CHECK(c.probes.size() == 1 && c.probes[0].name == "center" && c.probes[0].z == 30);
	CHECK(c.paths.output_dir == "./OUTPUT");

	// overrides: --set and -o / -i precedence
	{
		c = fdtd::config::load_sim_config(fdtd::config::build_config(BASIC, "x.ini", {"grid.dx=1e-9", "spectrum.points=51", "probe.center.x=11", "paths.output_dir=/ignored"},
												   "/data", "/scratch/run"));
		fdtd::config::validate_sim_config(c);
	}
	CHECK(c.grid.dx == 1e-9);
	CHECK(c.spectrum.points == 51);
	CHECK(c.probes[0].x == 11);
	CHECK(c.paths.output_dir == "/scratch/run"); // -o beats --set
	CHECK(c.paths.input_dir == "/data");

	// monitor boxes and superposed plane waves are INI sections
	c = load(std::string(BASIC) + "[monitor.gap]\ntype = time\ncenter_x = 5\ncenter_y = 6\ncenter_z = 7\n"
								  "size_x = 2\nsize_y = 2\nsize_z = 2\nstep = 4\n");
	CHECK(c.monitors.size() == 1 && c.monitors[0].name == "gap" && c.monitors[0].type == "time" && c.monitors[0].step == 4);
	CHECK(throws<fdtd::config::validation_error>([] { load(std::string(BASIC) + "[monitor.gap]\ncenter_x = 5\n"); }, "missing key:  [monitor.gap] size_x"));
	const std::string wave = "[wave.a]\ndirection_x = 0\ndirection_y = 1\ndirection_z = 0\ne_z = 1\nh_x = 1\n";
	CHECK(throws<fdtd::config::validation_error>([&] { load(std::string(BASIC) + wave); }, "type = wave_superposition"));
	c = load(std::string(BASIC) + wave, {"--set", "source.type=wave_superposition"});
	CHECK(c.waves.size() == 1 && c.waves[0].e_z == 1 && c.wave_superposition.dispersion_correction);
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "source.type=wave_superposition"}); }, "at least one [wave"));

	c = load(std::string(BASIC) + "[dipole.a]\nx = 3\ny = 4\nz = 5\n", {"--set", "source.type=dipole"});
	CHECK(c.dipoles.size() == 1 && c.dipoles[0].z == 5 && c.dipoles[0].e_z == 1);
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "source.type=dipole"}); }, "at least one [dipole"));

	// a key in the section of an unused type or signal is an error
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "signal.raised_cosine.max_frequency=1e15"}); }, "[signal.raised_cosine] has keys set"));
	CHECK(throws<fdtd::config::validation_error>([] { load(std::string(BASIC) + "[dipole.a]\nx = 3\ny = 4\nz = 5\n", {"--set", "source.type=dipole", "--set", "plane_wave.polarization_axis=x"}); },
										"[plane_wave] has keys set"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "plane_wave.polarization_axis=rhcp", "--set", "plane_wave.polarization_angle=10"}); }, "has no meaning"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "plane_wave.propagation_angle=20"}); }, "propagation_angle must be 0"));
	CHECK(throws<fdtd::config::syntax_error>([] { load(BASIC, {"--set", "plane_wave.propagation_axis=x"}); }, "expected one of"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "source.polarization=0"}); }, "unknown key"));
	c = load(BASIC, {"--set", "plane_wave.polarization_axis=x", "--set", "plane_wave.polarization_angle=30"});
	CHECK(c.plane_wave.polarization_angle == 30);
	c = fdtd::config::load_sim_config(fdtd::config::build_config(fdtd::config::to_ini(load(BASIC)), "resolved.ini", {})); // the resolved dump reloads
	fdtd::config::validate_sim_config(c);
	CHECK(c.signal.gaussian.center_wavelength == 800e-9);

	// shapes with a section of their own
	c = load(BASIC, {"--set", "geometry.shape=dipole", "--set", "dipole.width=2e-8", "--set", "dipole.top_end=4e-8", "--set", "dipole.bottom_end=3e-8"});
	{
		const fdtd::config::GeometrySection g = fdtd::config::resolved_geometry(c); // the builders read the pool
		CHECK(g.width == 2e-8 && g.length == 4e-8 && g.length_2 == 3e-8 && c.geometry.width == 0);
	}
	c = load(BASIC, {"--set", "geometry.shape=round_dipole", "--set", "round_dipole.rotation=30"});
	CHECK(fdtd::config::resolved_geometry(c).rotation == 120); // the builder counts from +x, the key from +z
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "geometry.shape=dipole", "--set", "geometry.width=2e-8"}); }, "unknown key")); // no generic dimensions
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "geometry.shape=dipole", "--set", "sphere.radius=2e-8"}); }, "[sphere] radius is set"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "geometry.shape=dipole", "--set", "dipole.radius=2e-8"}); }, "unknown key"));
	c = load(BASIC, {"--set", "geometry.shape=round_hybrid", "--set", "round_hybrid.width=2e-8", "--set", "round_hybrid.gap=1e-8", "--set", "round_hybrid.arm_length=3e-8"});
	{
		const fdtd::config::GeometrySection g = fdtd::config::resolved_geometry(c); // the builder's arm: (length_2 + arm_width + width - gap) / 2
		CHECK(std::fabs(0.5 * (g.length_2 + g.arm_width + g.width - g.gap) - 3e-8) < 1e-20);
	}
	CHECK(throws<fdtd::config::validation_error>([] { fdtd::config::shape_code("bullseye"); }, "unknown [geometry] shape")); // removed with patch and z_shape
	c = load(BASIC, {"--set", "geometry.shape=ellipse", "--set", "ellipse.radius_x=3e-8", "--set", "ellipse.radius_z=1e-8", "--set", "ellipse.rotation=25"});
	CHECK(fdtd::config::resolved_geometry(c).width == 3e-8 && fdtd::config::resolved_geometry(c).length == 1e-8 && fdtd::config::resolved_geometry(c).rotation == 25);
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "geometry.shape=cylinder_shell", "--set", "cylinder_shell.radius_core=2e-8", "--set", "cylinder_shell.radius_shell=1e-8"}); },
										"radius_shell > radius_core"));
	CHECK(fdtd::config::shape_code("cylinder_shell") == 115 && fdtd::config::find_shape_section("cylinder_shell") != nullptr);
	for (const char* shape : {"l_antenna", "lateral_dipole", "round_butterfly", "round_tripole", "x_shape"})
		CHECK(fdtd::config::find_shape_section(shape) != nullptr);
	c = load(BASIC, {"--set", "geometry.shape=round_tripole", "--set", "round_tripole.top_gap=2e-8", "--set", "round_tripole.bottom_gap=1e-8", "--set", "round_tripole.rotation=10"});
	{
		const fdtd::config::GeometrySection g = fdtd::config::resolved_geometry(c);
		CHECK(g.arm_length == 2e-8 && g.gap == 1e-8 && g.rotation == 100);
	}
	c = load(BASIC, {"--set", "geometry.shape=round_butterfly", "--set", "round_butterfly.length_z=3e-8"});
	CHECK(fdtd::config::resolved_geometry(c).length_2 == 3e-8 && fdtd::config::resolved_geometry(c).rotation == 45); // bars along x and z at 0
	c = load(BASIC, {"--set", "geometry.shape=x_shape", "--set", "x_shape.angle=60", "--set", "x_shape.rotated_width=1e-8"});
	CHECK(fdtd::config::resolved_geometry(c).rotation == 60 && fdtd::config::resolved_geometry(c).gap == 1e-8);

	// mode source and blocks
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "source.type=mode"}); }, "[mode] efield and hfield"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "source.type=mode", "--set", "mode.efield=e.h5", "--set", "mode.hfield=h.h5"}); },
										"effective_index"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "mode.position=1e-6"}); }, "[mode] has keys set"));
	CHECK(throws<fdtd::config::validation_error>([] {
		load(BASIC, {"--set", "source.type=mode", "--set", "mode.efield=e.h5", "--set", "mode.hfield=h.h5", "--set", "mode.effective_index=1.7",
					 "--set", "mode.propagation_axis=x", "--set", "mode.mpb_x_axis=x", "--set", "mode.mpb_y_axis=y"});
	}, "across the propagation axis"));
	CHECK(throws<fdtd::config::syntax_error>([] { load(BASIC, {"--set", "source.signal=mode"}); }, "expected one of")); // the mode is a source type now
	c = load(std::string(BASIC) + "[block.b]\nmaterial = 2\ny_max = 0\n[block.a]\nmaterial = 3\ncount = 4\npitch_x = 1e-6\n");
	CHECK(c.blocks.size() == 2 && c.blocks[0].name == "b" && c.blocks[0].y_max == 0 && c.blocks[0].x_min == -1 && c.blocks[1].count == 4); // file order
	CHECK(throws<fdtd::config::validation_error>([] { load(std::string(BASIC) + "[block.b]\nmaterial = 2\ny_min = 1e-6\ny_max = 0\n"); }, "min < max"));

	// typos and missing keys are errors, not silently ignored
	CHECK(throws<fdtd::config::validation_error>([] { load(std::string(BASIC) + "[grid2]\nfoo = 1\n"); }, "unknown key:  [grid2] foo"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "grid.dX_typo=3"}); }, "dx_typo"));
	CHECK(throws<fdtd::config::validation_error>([] { load("[grid]\nnx=1\n"); }, "missing key:  [grid] dx"));
	CHECK(throws<fdtd::config::syntax_error>([] { load(std::string(BASIC) + "[run]\ntrace_interval = ten\n"); }, "trace_interval"));
	CHECK(throws<fdtd::config::syntax_error>([] { load(BASIC, {"--set", "source.signal=triangle"}); }, "expected one of"));
	CHECK(throws<fdtd::config::syntax_error>([] { load("nx = 1\n"); }, "outside of a [section]"));

	// validation
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "geometry.first_medium=7"}); }, "[materials] n7 is required"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "parallel.chunks_x=0 0"}); }, "strictly increasing"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "parallel.chunks_x=0"}); }, "exactly 2 entries"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "nonlinear.shg=antenna", "--set", "spectrum.harmonic_band=second", "--set", "spectrum.spacing=frequency"}); }, "chi_structure"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "multipoles.sphere_radius=4"}); }, "set together"));
	c = load(BASIC, {"--set", "parallel.chunks_x=0 120"});

	// nonlinear tensors: [chi2|chi3.<name>] sections and presets
	{
		const char* shg[] = {"--set", "nonlinear.shg=antenna", "--set", "spectrum.harmonic_band=second", "--set", "spectrum.spacing=frequency"};
		std::vector<const char*> on(std::begin(shg), std::end(shg));
		auto with = [&](const char* k) { auto v = on; v.push_back("--set"); v.push_back(k); return v; };
		const std::string mine = std::string(BASIC) + "[chi2.Mine]\nd22 = 2.1e-12\nd31 = 0.15e-12\n[chi3.Au]\nxxxx = 1\n";
		auto t = load(mine, with("nonlinear.chi_structure=mine"));
		const fdtd::config::ChiTensor* m = fdtd::config::find_chi_tensor(2, "MINE", t.chi_tensors);
		CHECK(m && m->values[7] == 2.1e-12 && m->values[12] == 0.15e-12 && m->values[0] == 0);
		CHECK(fdtd::config::find_chi_tensor(3, "au", t.chi_tensors)->values[1] == 0);                   // a section replaces the preset
		CHECK(fdtd::config::find_chi_tensor(3, "au", {})->values[1] == 2.533e-19);
		CHECK(fdtd::config::find_chi_tensor(2, "BBO", {})->values[7] == 2.1e-12);
		CHECK(throws<fdtd::config::validation_error>([&] { load(mine, with("nonlinear.chi_structure=other")); }, "no [chi2.other] section"));
		CHECK(throws<fdtd::config::validation_error>([&] { load(mine, with("nonlinear.chi_structure=Au")); }, "chi2 presets: Zero, BBO"));
		CHECK(throws<fdtd::config::syntax_error>([&] { load(std::string(BASIC) + "[chi2.x]\nxxxx = 1\n"); }, "d11 ... d36"));
		CHECK(fdtd::config::to_ini(load(fdtd::config::to_ini(t))) == fdtd::config::to_ini(t));
	}
	CHECK(fdtd::config::parse_int_list("d", "k", c.parallel.chunks_x) == std::vector<int>({0, 120}));


	// spectrum: the solver gets exactly the requested frequencies
	{
		const double c0 = 299792458.0;
		const double pi2 = 6.283185307179586;
		auto close = [](double a, double b) { return std::fabs(a - b) <= 1e-12 * std::fabs(b); };
		auto w = fdtd::config::spectrum_omegas(load(BASIC, {}).spectrum);                  // 41 points uniform in wavelength
		auto s0 = load(BASIC, {}).spectrum;
		CHECK(w.size() == 41 && close(w.front(), pi2 * c0 / s0.lambda_max) && close(w.back(), pi2 * c0 / s0.lambda_min));
		w = fdtd::config::spectrum_omegas(load(BASIC, {"--set", "spectrum.points=1", "--set", "spectrum.lambda_max=500e-9"}).spectrum);
		CHECK(w.size() == 1 && close(w[0], pi2 * c0 / s0.lambda_min));
		const char* shg[] = {"--set", "nonlinear.shg=bulk", "--set", "nonlinear.chi_structure=BBO", "--set", "nonlinear.chi_first=zero",
							 "--set", "nonlinear.chi_middle=Zero", "--set", "nonlinear.chi_last=Zero"};
		std::vector<const char*> base(std::begin(shg), std::end(shg));
		auto with = [&](std::vector<const char*> extra) { auto v = base; v.insert(v.end(), extra.begin(), extra.end()); return v; };
		w = fdtd::config::spectrum_omegas(load(BASIC, with({"--set", "spectrum.harmonic_band=second", "--set", "spectrum.points=1", "--set", "spectrum.lambda_max=500e-9"})).spectrum); // pair
		CHECK(w.size() == 2 && close(w[1], 2 * w[0]));
		w = fdtd::config::spectrum_omegas(load(BASIC, with({"--set", "spectrum.harmonic_band=second", "--set", "spectrum.spacing=frequency"})).spectrum);
		CHECK(close(w[1] - w[0], w[2] - w[1]) && close(w[41], 2 * w[0]) && close(w[81], 2 * w[40])); // band: uniform in frequency
		CHECK(fdtd::config::spectrum_omegas(load(BASIC, with({"--set", "spectrum.harmonic_band=second", "--set", "spectrum.spacing=frequency"})).spectrum).size() == 82);
		// no mystery mutations: a harmonic without an explicit band is rejected, not silently doubled
		CHECK(throws<fdtd::config::validation_error>([&] { load(BASIC, base); }, "harmonic_band = second"));
		CHECK(throws<fdtd::config::validation_error>([&] { load(BASIC, {"--set", "spectrum.harmonic_band=third"}); }, "requires [nonlinear] thg"));
		CHECK(throws<fdtd::config::validation_error>([&] { load(BASIC, with({"--set", "spectrum.harmonic_band=second"})); }, "spacing = frequency"));
		CHECK(throws<fdtd::config::validation_error>([&] { load(BASIC, {"--set", "spectrum.points=1"}); }, "lambda_min = lambda_max"));
	}

	// resolved dump round-trips
	c = load(BASIC, {"--set", "grid.dx=1.234567890123e-9"});
	auto again = load(fdtd::config::to_ini(c));
	CHECK(fdtd::config::to_ini(again) == fdtd::config::to_ini(c));
	CHECK(again.grid.dx == 1.234567890123e-9);
	CHECK(again.probes.size() == 1 && again.materials.refractive_index.size() == 1);

	// template parses, all keys are known
	CHECK(throws<fdtd::config::syntax_error>([] { fdtd::config::Config t; t.load_string(fdtd::config::template_ini() + "\n[grid]\n"); }, "duplicate section"));
	{
		// the template with the required keys filled in must load cleanly
		std::string t = fdtd::config::template_ini();
		for (const char* k : {"nx", "ny", "nz"}) { auto p = t.find(std::string("\n") + k + " = 0"); t.replace(p, 6 + std::string(k).size() - 2, std::string("\n") + k + " = 10"); }
		for (const char* k : {"dx", "dy", "dz"}) { auto p = t.find(std::string("\n") + k + " = 0"); t.replace(p, 6 + std::string(k).size() - 2, std::string("\n") + k + " = 1e-9"); }
		{ auto p = t.find("\npoints = 1"); t.replace(p, 11, "\npoints = 5"); }
		bool ok = true;
		try { load(t); } catch (const std::exception& e) { ok = false; std::printf("  template: %s\n", e.what()); }
		CHECK(ok);
	}

	// overrides
	CHECK(throws<std::invalid_argument>([&] { fdtd::config::build_config(BASIC, "x.ini", {"nodot=1"}); }, "SECTION.KEY=VALUE"));
	CHECK(fdtd::config::load_sim_config(fdtd::config::build_config(BASIC, "x.ini", {}, "in", "out")).paths.output_dir == "out");

	// named keys that replaced numeric codes
	CHECK(fdtd::config::shape_code("none") == 0 && fdtd::config::shape_code("sphere") == 1 && fdtd::config::shape_code("Sphere_Dimer") == 67 &&
		  fdtd::config::shape_code("rounded_dipole_offset_planes") == 114);
	CHECK(throws<fdtd::config::validation_error>([] { fdtd::config::shape_code("5"); }, "unknown [geometry] shape"));
	CHECK(throws<fdtd::config::syntax_error>([] { load(BASIC, {"--set", "geometry.shape=spherical"}); }, "expected one of"));
	{
		fdtd::config::OutputSection o;
		CHECK(fdtd::config::dft_plane_mask_bits(o) == 0);
		o.dft_yz = true;
		CHECK(fdtd::config::dft_plane_mask_bits(o) == 1); // legacy DFT Plot 1
		o.dft_xz = o.dft_xy = true;
		CHECK(fdtd::config::dft_plane_mask_bits(o) == 7); // legacy 4: all planes
		o.dft_magnetic = true;
		CHECK(fdtd::config::dft_plane_mask_bits(o) == 15); // legacy 5: all planes + H
	}
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "output.dft_magnetic=true"}); }, "dft_magnetic needs"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "hydrodynamic.pressure=false"}); }, "[hydrodynamic] has keys set"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "two_temperature.save_frames=true"}); }, "[two_temperature] has keys set"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "pump_probe.delay=1e-13", "--set", "temperature.model=two_temperature"}); }, "[pump_probe] has keys set"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "electrons.diffusion=1"}); }, "need model = nonlocal"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "anapole.patch_edge=1e-7"}); }, "[anapole] has keys set"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "electrons.model=gnor", "--set", "temperature.model=two_temperature"}); },
										"models combine"));
	CHECK(throws<fdtd::config::syntax_error>([] { load(BASIC, {"--set", "electrons.model=drude"}); }, "expected one of"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "model.electrons=gnor"}); }, "unknown key"));
	CHECK(load(BASIC, {"--set", "temperature.model=two_temperature", "--set", "two_temperature.save_frames=true"}).two_temperature.save_frames);

	// rank split planner
	CHECK((fdtd::config::rank_starts("x", 52, 4, "", {16, 36}) == std::vector<int>{0, 13, 26, 39}));
	{
		std::string note;
		auto s6 = fdtd::config::rank_starts("x", 52, 6, "", {16, 36}, &note); // equal split puts a boundary on 36
		CHECK((s6 == std::vector<int>{0, 9, 18, 27, 37, 44}) && note.find("moved to 37") != std::string::npos);
		auto s7 = fdtd::config::rank_starts("y", 52, 7, "", {16, 36}); // ... and on 16
		CHECK(std::find(s7.begin(), s7.end(), 16) == s7.end() && s7.size() == 7);
	}
	CHECK(throws<fdtd::config::validation_error>([] { fdtd::config::rank_starts("x", 52, 2, "0 16", {16, 36}); }, "lies on the TF/SF plane"));
	CHECK(throws<fdtd::config::validation_error>([] { fdtd::config::rank_starts("x", 52, 27, "", {16, 36}); }, "at least 2"));
	CHECK(throws<fdtd::config::validation_error>([] { fdtd::config::rank_starts("x", 52, 26, "", {16, 36}); }, "no valid split"));
	CHECK(fdtd::config::rank_starts("x", 52, 1, "", {16, 36}) == std::vector<int>{0});

	std::printf(failures ? "%d FAILURE(S)\n" : "all config/CLI tests passed\n", failures);
	return failures ? 1 : 0;
}
