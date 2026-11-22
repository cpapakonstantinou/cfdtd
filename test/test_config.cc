#include <algorithm>
// Unit tests for fdtd::config::Config, fdtd::config::SimConfig and the CLI layer, no MPI needed.
#include "config/errors.h"
#include "config/sim_config.h"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
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

[geometry]
material = Glass

[material.glass]
model = dielectric
n = 1.45

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
	CHECK(c.definitions.size() == 2 && c.definitions[0].name == "vacuum" && c.definitions[1].name == "glass" && c.definitions[1].n == 1.45);
	CHECK(fdtd::config::material_slot(c, "Glass") == 1 && fdtd::config::material_slot(c, "vacuum") == 0);
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
	c = load(std::string(BASIC) + "[block.b]\nmaterial = glass\ny_max = 0\n[block.a]\nmaterial = vacuum\ncount = 4\npitch_x = 1e-6\n");
	CHECK(c.blocks.size() == 2 && c.blocks[0].name == "b" && c.blocks[0].y_max == 0 && c.blocks[0].x_min == -1 && c.blocks[1].count == 4); // file order
	CHECK(throws<fdtd::config::validation_error>([] { load(std::string(BASIC) + "[block.b]\nmaterial = glass\ny_min = 1e-6\ny_max = 0\n"); }, "min < max"));

	// typos and missing keys are errors, not silently ignored
	CHECK(throws<fdtd::config::validation_error>([] { load(std::string(BASIC) + "[grid2]\nfoo = 1\n"); }, "unknown key:  [grid2] foo"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "grid.dX_typo=3"}); }, "dx_typo"));
	CHECK(throws<fdtd::config::validation_error>([] { load("[grid]\nnx=1\n"); }, "missing key:  [grid] dx"));
	CHECK(throws<fdtd::config::syntax_error>([] { load(std::string(BASIC) + "[run]\ntrace_interval = ten\n"); }, "trace_interval"));
	CHECK(throws<fdtd::config::syntax_error>([] { load(BASIC, {"--set", "source.signal=triangle"}); }, "expected one of"));
	CHECK(throws<fdtd::config::syntax_error>([] { load("nx = 1\n"); }, "outside of a [section]"));

	// validation
	// materials: definitions, extends, poles, errors
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "geometry.first_medium=water"}); }, "material 'water' is not defined"));
	CHECK(throws<fdtd::config::syntax_error>([] { load(std::string(BASIC) + "[material.bad]\nmodel = dielectric\nn = 1\nplasma_frequency = 1\n"); }, "not a key of model dielectric"));
	CHECK(throws<fdtd::config::syntax_error>([] { load(std::string(BASIC) + "[material.bad]\nmodel = plasma\n"); }, "not a model"));
	CHECK(throws<fdtd::config::validation_error>([] { load(std::string(BASIC) + "[material.vacuum]\nmodel = dielectric\nn = 2\n"); }, "built in"));
	{
		const std::string metal = std::string(BASIC) +
			"[material.metal]\nmodel = drude_cp\neps_inf = 1.5\nplasma_frequency = 1e16\nplasma_damping = 1e14\n"
			"[material.metal.pole.1]\namplitude = 1\nfrequency = 4e15\ndamping = 5e14\nphase = -1\n"
			"[material.metal.pole.2]\namplitude = 2\nfrequency = 5e15\ndamping = 6e14\nphase = 1\n"
			"[material.wet]\nextends = metal\nplasma_damping = 2e14\n[material.wet.pole.2]\nphase = 0.5\n";
		auto m = load(metal, {"--set", "geometry.material=wet", "--set", "geometry.last_medium=metal"});
		const auto& w = m.definitions[fdtd::config::material_slot(m, "wet")];
		CHECK(w.kind == fdtd::materials::Kind::drude_cp && w.eps_inf == 1.5 && w.plasma_damping == 2e14 && w.poles.size() == 2);
		CHECK(w.poles[1].phase == 0.5 && w.poles[1].amplitude == 2 && w.poles[0].phase == -1);
		CHECK(fdtd::config::material_slot(m, "metal") == 2 && m.definitions.size() == 3);
		CHECK(fdtd::config::to_ini(load(fdtd::config::to_ini(m))) == fdtd::config::to_ini(m));
		CHECK(throws<fdtd::config::validation_error>([&] { load(metal + "[material.metal.pole.4]\namplitude = 1\n", {"--set", "geometry.material=metal"}); }, "without a gap"));
		CHECK(throws<fdtd::config::validation_error>([&] { load(metal + "[material.loop]\nextends = loop\n", {"--set", "geometry.material=loop"}); }, "circle"));
	}
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "parallel.chunks_x=0 0"}); }, "strictly increasing"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "parallel.chunks_x=0"}); }, "exactly 2 entries"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "nonlinear.shg=true", "--set", "spectrum.harmonic_band=second", "--set", "spectrum.spacing=frequency"}); }, "a material with a chi2 tensor"));
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "multipoles.sphere_radius=4"}); }, "set together"));
	c = load(BASIC, {"--set", "parallel.chunks_x=0 120"});

	// nonlinear tensors: [chi2|chi3.<name>] sections named by a material
	{
		const char* shg[] = {"--set", "nonlinear.shg=true", "--set", "spectrum.harmonic_band=second", "--set", "spectrum.spacing=frequency"};
		std::vector<const char*> on(std::begin(shg), std::end(shg));
		const std::string mine = std::string(BASIC) + "[chi2.Mine]\nd22 = 2.1e-12\nd31 = 0.15e-12\n[material.crystal]\nextends = glass\nchi2 = Mine\n";
		auto t = load(mine, {"--set", "geometry.material=crystal", "--set", "nonlinear.shg=true", "--set", "spectrum.harmonic_band=second", "--set", "spectrum.spacing=frequency"});
		const auto& d = t.definitions[1];
		CHECK(d.name == "crystal" && d.chi2 == "mine" && d.chi2_values[7] == 2.1e-12 && d.chi2_values[12] == 0.15e-12 && d.chi2_values[0] == 0 && d.n == 1.45);
		CHECK(throws<fdtd::config::validation_error>([&] { load(mine + "[material.other]\nextends = glass\nchi2 = none_such\n", {"--set", "geometry.material=other"}); }, "no [chi2.none_such] section"));
		CHECK(throws<fdtd::config::syntax_error>([&] { load(std::string(BASIC) + "[chi2.x]\nxxxx = 1\n[material.x]\nextends = glass\nchi2 = x\n", {"--set", "geometry.material=x"}); }, "d11 ... d36"));
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
		const char* shg[] = {"--set", "nonlinear.shg=true", "--set", "geometry.material=crystal", "--set", "material.crystal.extends=glass",
							 "--set", "material.crystal.chi2=c", "--set", "chi2.c.d11=1e-12"};
		std::vector<const char*> base(std::begin(shg), std::end(shg));
		auto with = [&](std::vector<const char*> extra) { auto v = base; v.insert(v.end(), extra.begin(), extra.end()); return v; };
		w = fdtd::config::spectrum_omegas(load(BASIC, with({"--set", "spectrum.harmonic_band=second", "--set", "spectrum.points=1", "--set", "spectrum.lambda_max=500e-9"})).spectrum); // pair
		CHECK(w.size() == 2 && close(w[1], 2 * w[0]));
		w = fdtd::config::spectrum_omegas(load(BASIC, with({"--set", "spectrum.harmonic_band=second", "--set", "spectrum.spacing=frequency"})).spectrum);
		CHECK(close(w[1] - w[0], w[2] - w[1]) && close(w[41], 2 * w[0]) && close(w[81], 2 * w[40])); // band: uniform in frequency
		CHECK(fdtd::config::spectrum_omegas(load(BASIC, with({"--set", "spectrum.harmonic_band=second", "--set", "spectrum.spacing=frequency"})).spectrum).size() == 82);
		// no mystery mutations: a harmonic without an explicit band is rejected, not silently doubled
		CHECK(throws<fdtd::config::validation_error>([&] { load(BASIC, base); }, "harmonic_band = second"));
		// a band without a nonlinearity is a linear run that also samples k w; a band other than the process is rejected
		CHECK(fdtd::config::spectrum_omegas(load(BASIC, {"--set", "spectrum.harmonic_band=third"}).spectrum).size() == 82);
		CHECK(throws<fdtd::config::validation_error>([&] { load(BASIC, with({"--set", "spectrum.harmonic_band=third"})); }, "does not match [nonlinear] shg"));
		w = fdtd::config::spectrum_omegas(load(BASIC, with({"--set", "spectrum.harmonic_band=second"})).spectrum); // band uniform in wavelength
		CHECK(w.size() == 82 && !close(w[1] - w[0], w[2] - w[1]) && close(w[41], 2 * w[0]) && close(w[81], 2 * w[40]));
		CHECK(throws<fdtd::config::validation_error>([&] { load(BASIC, {"--set", "spectrum.points=1"}); }, "lambda_min = lambda_max"));

		// Gauss-Legendre: nodes inside the band, symmetric in frequency, odd counts at its middle; exact for
		// polynomials of degree 2n - 1; the trapezoid weights sum to the band too
		for (int n : {2, 5, 8})
		{
			const char* gauss[] = {"--set", "spectrum.spacing=gauss", "--set", "spectrum.lambda_min=400e-9", "--set", "spectrum.lambda_max=800e-9"};
			std::vector<const char*> args(std::begin(gauss), std::end(gauss));
			const std::string points = "spectrum.points=" + std::to_string(n);
			args.push_back("--set");
			args.push_back(points.c_str());
			const auto s = load(BASIC, args).spectrum;
			const double two_pi = 6.283185307179586;
			const std::vector<double> w = fdtd::config::spectrum_omegas(s);
			const std::vector<double> q = fdtd::config::spectrum_weights(s);
			const double a = 299792458.0 / 800e-9;
			const double b = 299792458.0 / 400e-9;
			double sum = 0;
			double moment = 0;
			bool inside = true;
			for (int m = 0; m < n; ++m)
			{
				const double f = w[m] / two_pi;
				const double u = (f - a) / (b - a);
				inside &= f > a && f < b && (m == 0 || w[m] > w[m - 1]);
				sum += q[m];
				moment += q[m] * std::pow(u, 2 * n - 1);
			}
			CHECK(inside && close(sum, b - a) && std::fabs(moment / (b - a) - 1.0 / (2 * n)) < 1e-12);
			if (n == 5)
				CHECK(close(w[2] / two_pi, 0.5 * (a + b)));
		}
		{
			const auto s = load(BASIC, {"--set", "spectrum.spacing=wavelength"}).spectrum;
			double sum = 0;
			for (double v : fdtd::config::spectrum_weights(s))
				sum += v;
			CHECK(close(sum, 299792458.0 / s.lambda_min - 299792458.0 / s.lambda_max));
		}
	}

	// the catalog: files outside the binary, named by the materials and tensors they define
	{
		const std::string dir = "/tmp/fdtd_test_catalog";
		std::filesystem::create_directories(dir + "/chi3");
		std::ofstream(dir + "/metal.ini") << "[material.metal]\nmodel = drude\neps_inf = 1\nplasma_frequency = 1e16\nplasma_damping = 1e14\nchi3 = kerr\n";
		std::ofstream(dir + "/chi3/kerr.ini") << "[chi3.kerr]\nxxxx = 1e-19\n";
		const std::string run = std::string(BASIC) + "[materials]\ncatalog = " + dir + "\n";
		const std::vector<std::string> set = {"geometry.last_medium=metal", "material.mine.extends=metal", "material.mine.plasma_damping=2e14", "geometry.first_medium=mine"};
		const std::string text = fdtd::config::with_catalog(run, "x.ini", set);
		CHECK(text.find("[material.metal]") != std::string::npos && text.find("[chi3.kerr]") != std::string::npos);
		auto m = fdtd::config::load_sim_config(fdtd::config::build_config(text, "x.ini", set));
		CHECK(m.definitions.size() == 4 && m.definitions[2].name == "mine" && m.definitions[2].plasma_damping == 2e14 && m.definitions[3].chi3_values[0] == 1e-19);
		CHECK(throws<std::invalid_argument>([&] { fdtd::config::with_catalog(run, "x.ini", {"geometry.material=unobtainium"}); }, "neither in the run file nor in the catalog"));
	}

	// resolved dump round-trips
	c = load(BASIC, {"--set", "grid.dx=1.234567890123e-9"});
	auto again = load(fdtd::config::to_ini(c));
	CHECK(fdtd::config::to_ini(again) == fdtd::config::to_ini(c));
	CHECK(again.grid.dx == 1.234567890123e-9);
	CHECK(again.probes.size() == 1 && again.definitions.size() == 2 && again.definitions[1] == c.definitions[1]);

	// template parses, all keys are known
	CHECK(throws<fdtd::config::syntax_error>([] { fdtd::io::Ini t; t.load_string(fdtd::config::template_ini() + "\n[grid]\n"); }, "duplicate section"));
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
	CHECK(throws<fdtd::config::validation_error>([] { load(BASIC, {"--set", "electrons.fermi_velocity=2"}); }, "needs model = nonlocal"));
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
