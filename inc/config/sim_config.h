#ifndef __SIM_CONFIG_H__
#define __SIM_CONFIG_H__
/** \file   config/sim_config.h
 *  \brief  Typed simulation configuration built on top of Config, INI.
 *
 *  Rules of the interface:
 *    - one flat scalar per key, snake_case, 1-2 self-documenting words
 *      [grid] dx / dy / dz, never "Space Step = [dx, dy, dz]";
 *    - choices are names, pulse = gaussian, never magic integers;
 *    - a key means exactly one thing: no key changes the meaning of another
 *      the solver takes its DFT frequencies straight from spectrum_omegas();
 *    - unknown keys are errors.
 *
 *  The schema is declared exactly once, in visit_sim_config(), simConfig.cc,
 *  which drives loading, unknown-key detection, the resolved dump and the
 *  `fdtd template` output. No dependency on MPI or the solver globals.
 **/
#include <map>
#include <string>
#include <vector>
#include "config/chi_tensors.h"
#include "config/ini.h"
#include "geometry/anapole.h"
#include "geometry/block.h"
#include "sources/mode_spec.h"
#include "monitors/inputs.h"

namespace fdtd::config
{

struct RunSection
{
	double wallclock_limit_hours = 0.0; ///< stop the time loop before this wall time [h]; 0 = never, was hardcoded 23.5 h
	int wallclock_check_step  = 1000; ///< step at which the wall-clock limit is extrapolated
	int trace_interval        = 100; ///< print the field at the animation point every N steps; 0 = off
	int checkpoint_count      = 10; ///< probes / trackers / Q are flushed this many times per run
};

struct PathsSection
{
	std::string input_dir  = ".";        ///< data files, mode profiles, electron temperature
	std::string output_dir = "./OUTPUT"; ///< root of all output
};

struct GridSection
{
	int nx = 0; ///< interior cells, PML is added on top
	int ny = 0; ///< interior cells, PML is added on top
	int nz = 0; ///< interior cells, PML is added on top
	double dx = 0; ///< cell size [m]
	double dy = 0; ///< cell size [m]
	double dz = 0; ///< cell size [m]
	int pml_x = 20; ///< CPML thickness per face [cells]
	int pml_y = 20; ///< CPML thickness per face [cells]
	int pml_z = 20; ///< CPML thickness per face [cells]
	int time_steps = 100000;
	double courant = 0.5;                ///< Courant number used to derive dt
	int origin_x = -1; ///< structure origin [cells]; -1 = centre
	int origin_y = -1; ///< structure origin [cells]; -1 = centre
	int origin_z = -1; ///< structure origin [cells]; -1 = centre
	bool periodic_x = false;
	bool periodic_z = false;
	bool staircase = false; ///< staircased instead of per-component material assignment
};

struct ParallelSection
{
	int ranks_x = 1; ///< MPI ranks per axis; product = MPI world size
	int ranks_y = 1; ///< MPI ranks per axis; product = MPI world size
	int ranks_z = 1; ///< MPI ranks per axis; product = MPI world size
	std::string chunks_x; ///< optional chunk start indices ("0 40 90"); empty = equal split
	std::string chunks_y; ///< optional chunk start indices ("0 40 90"); empty = equal split
	std::string chunks_z; ///< optional chunk start indices ("0 40 90"); empty = equal split
};

struct BoxesSection
{
	int tfsf_x = 11; ///< TF/SF box margin from the PML [cells]
	int tfsf_y = 11; ///< TF/SF box margin from the PML [cells]
	int tfsf_z = 11; ///< TF/SF box margin from the PML [cells]
	int scattering_x = -1; ///< scattering box margin; -1 = tfsf - 5
	int scattering_y = -1; ///< scattering box margin; -1 = tfsf - 5
	int scattering_z = -1; ///< scattering box margin; -1 = tfsf - 5
	int absorption_x = -1; ///< absorption box margin; -1 = tfsf + 5
	int absorption_y = -1; ///< absorption box margin; -1 = tfsf + 5
	int absorption_z = -1; ///< absorption box margin; -1 = tfsf + 5
};

/** [source]: what is injected, type, and its time signal. The keys of one source type are in the
 *  section named after it, [plane_wave], [gaussian_beam], [wave_superposition], [dipole.<name>], and the
 *  keys of one signal in [signal.<name>]. Setting a key in the section of a type or signal that the run
 *  does not use is an error. */
struct SourceSection
{
	std::string type = "plane_wave";       ///< plane_wave | dipole | gaussian_beam | wave_superposition | mode
	std::string signal = "raised_cosine";  ///< see pulse_code()
	double peak_field = 1.0;               ///< [V/m]
	bool allow_truncated_signal = false; ///< run even when the time steps do not cover the signal
	int aux_steps = -1; ///< 1D auxiliary grid (plane_wave, wave_superposition): steps; -1 = [grid] time_steps
	int aux_length_factor = 1; ///< 1D auxiliary grid: length as a multiple of the 3D y extent
	bool operator==(const SourceSection&) const = default;
};

/** [plane_wave]: direction and polarisation of the plane wave.
 *  The wave travels along +propagation_axis, only y for now. polarization_axis names the transverse
 *  axis of E (x or z) and polarization_angle turns E from it toward the other transverse axis;
 *  rhcp / lhcp give circular polarisation, IEEE convention: seen along the propagation direction. */
struct PlaneWaveSection
{
	std::string propagation_axis = "y";    ///< y
	double propagation_angle = 0.0;        ///< [deg] 0, oblique incidence: wave_superposition
	std::string polarization_axis = "z";   ///< x | z | rhcp | lhcp
	int polarization_angle = 0; ///< [deg] from polarization_axis toward the other transverse axis
	bool operator==(const PlaneWaveSection&) const = default;
};

/// [wave_superposition]: options of the waves listed in the [wave.<name>] sections
struct WaveSuperpositionSection
{
	bool derivative = false;               ///< inject the time derivative of the signal
	bool dispersion_correction = true;     ///< correct each wave for grid dispersion
	bool layered_simulation = false;       ///< run the 1D simulations of the waves even when the three layers are one medium
	bool operator==(const WaveSuperpositionSection&) const = default;
};

/// [signal.raised_cosine], [signal.raised_cosine_derivative], [signal.rectangular]
struct BandLimitedSignal
{
	double max_frequency = 2.0e15;         ///< [Hz] bandwidth
	bool operator==(const BandLimitedSignal&) const = default;
};

/// [signal.gaussian]: modulated Gaussian
struct GaussianSignal
{
	double center_wavelength = 1000e-9;    ///< [m] carrier
	double width = 5e-15;                  ///< [s]
	double chirp = 0.0;
	double chirp_quadratic = 0.0;
	bool operator==(const GaussianSignal&) const = default;
};

/// [signal.gaussian_alt]
struct GaussianAltSignal
{
	double center_wavelength = 1000e-9;    ///< [m] carrier
	double width = 5e-15;                  ///< [s]
	bool operator==(const GaussianAltSignal&) const = default;
};

/// [signal.continuous_wave]: the carrier is the centre of [spectrum]
struct ContinuousWaveSignal
{
	double settle_periods = 50.0;          ///< the DFT starts after this many periods
	bool operator==(const ContinuousWaveSignal&) const = default;
};

/// [signal.cw_pulse]: the carrier is the centre of [spectrum]
struct CwPulseSignal
{
	int periods = 30;                      ///< signal length [periods]
	bool operator==(const CwPulseSignal&) const = default;
};

/// the [signal.<name>] sections, constant and sinc have no keys
struct SignalSections
{
	BandLimitedSignal raised_cosine;
	BandLimitedSignal raised_cosine_derivative;
	BandLimitedSignal rectangular;
	GaussianSignal gaussian;
	GaussianAltSignal gaussian_alt;
	ContinuousWaveSignal continuous_wave;
	CwPulseSignal cw_pulse;
};

using ModeSection = fdtd::sources::ModeSpec; ///< [mode], sources/mode_spec.h

struct GaussianBeamSection
{
	int theta_steps = 0;
	int phi_steps = 0;
	double theta_max = 0; ///< [deg]
	double phi_max = 0; ///< [deg]
	double waist_factor = 0;
	double focal_factor = 0;
	bool operator==(const GaussianBeamSection&) const = default;
};

/** The DFT frequencies. lambda_min..lambda_max with `points` samples means
 *  exactly that, in every simulation. A harmonic band is an *additional*
 *  band at 2x / 3x the frequencies of the fundamental band, switched on
 *  explicitly; it never alters lambda_min, lambda_max or points. */
struct SpectrumSection
{
	double lambda_min = 400e-9;     ///< [m]
	double lambda_max = 800e-9;     ///< [m]
	int points = 1; ///< samples in [lambda_min, lambda_max]; 1 = single wavelength (lambda_min == lambda_max)
	std::string spacing = "wavelength"; ///< wavelength | frequency, uniform step in lambda or in omega
	std::string harmonic_band = "none";  ///< none | second | third : extra band at 2w / 3w
};

struct OutputSection
{
	bool dft_xy = false;           ///< frequency-domain field on the XY plot plane
	bool dft_xz = false;           ///< ... on the XZ plot plane
	bool dft_yz = false;           ///< ... on the YZ plot plane
	bool dft_magnetic = false;     ///< also write H on the enabled planes
	bool volume_dft = false;       ///< full-volume E/H DFT
	int frames_xy = 0; ///< time-domain frame period [steps]; 0 = off
	int frames_xz = 0; ///< time-domain frame period [steps]; 0 = off
	int frames_yz = 0; ///< time-domain frame period [steps]; 0 = off
	int frames_3d = 0; ///< time-domain frame period [steps]; 0 = off
	bool phase_normalize = true;   ///< normalize the DFT phase to the source
	bool sigma_absorption = false; ///< absorption from the sigma|E|^2 volume integral
	bool xdmf = true;              ///< write a .xmf next to every HDF5 output, ParaView / VisIt
};

struct MultipolesSection
{
	bool exact = false;            ///< exact multipole decomposition, Alaee et al. 2018
	bool lwa = false;              ///< long-wavelength multipoles, the kr -> 0 limit, from the same volume current
	int sphere_radius = 0; ///< legacy sphere-sampled multipoles: radius [cells]; 0 = off
	int sphere_theta = 0; ///< sampling points on that sphere
	int sphere_phi = 0; ///< sampling points on that sphere
};

struct GeometrySection
{
	std::string shape = "none";    ///< structure, see shape_code()
	int material = 0; ///< material of the structure
	int first_medium = 0;
	int middle_medium = 0;
	int last_medium = 0;
	double middle_thickness = 0.0; ///< [m]
	double shift_y = 0.0;          ///< structure offset along y [m]
	// the pool of dimensions of the shapes without a section of their own (shape_sections())
	// Not INI keys: the dimensions the builders read, filled from the shape's own section (resolved_geometry()).
	double width = 0; ///< [m]
	double length = 0; ///< [m]
	double thickness = 0; ///< [m]
	double length_2 = 0; ///< [m]
	double gap = 0; ///< [m]
	double arm_length = 0; ///< [m]
	double arm_width = 0; ///< [m]
	double rotation = 0;           ///< [deg]
	double corner_radius = 0;      ///< [m]
	int shell_material = 2;
};

struct MaterialsSection
{
	std::map<int, double> refractive_index; ///< n<i> = value, i = 0..99, dispersionless media
	double carrier_density = 0;
	double carrier_density_back = 0;
	double stretch_x = 0;           ///< graphene, material 102: scales the plasma frequency and the critical-point amplitudes
};

struct CustomMaterialSection
{
	double eps_inf = 0;
	double plasma_frequency = 0;
	double plasma_damping = 0;
	double cp_frequency = 0;
	double cp_frequency_2 = 0;
	double cp_damping = 0;
	double cp_damping_2 = 0;
	double cp_amplitude = 0;
	double cp_amplitude_2 = 0;
	double cp_phase = 0;
	double cp_phase_2 = 0;
	double fermi_velocity = 0;
	double electron_density = 0;
	bool operator==(const CustomMaterialSection&) const = default;
};

/**
 * @brief [anapole]: plasmonic anapole nanoantenna, shape = anapole
 *
 * A patch of edge `patch_edge` with rounded corners, minus a dumbbell aperture: two lobes joined
 * by a rectangular gap. u runs along the lobe axis, x at rotation 0, v across it (z).
 *
 * @image html anapole.png "[anapole] keys on the mesh, docs/geometry/Shapes.md"
 *
 * The lobe outline is v = +-patch_edge lobe_height s^lobe_exponent sqrt(1 - s^2), s = (|u| / patch_edge - lobe_offset) / lobe_length,
 * 0 < s <= 1: it starts at the gap face |u| = gap_length / 2 and ends at the tip |u| = patch_edge (lobe_offset + lobe_length).
 * Where the lobe is wider than the gap at the gap face, a metal tooth separates the two. The drawing is extruded over `thickness` along y.
 *
 * The lobe coefficients are in units of patch_edge and were fitted at patch_edge = 200 nm; the gap is absolute.
 * aperture_material fills the lobes and the gap, gap_material the gap alone.
 */
using AnapoleSection = geom::AnapoleSpec; ///< [anapole], geometry/anapole.h

struct NonlinearSection
{
	std::string shg = "none";       ///< none | antenna | bulk : second-harmonic generation in the structure or in the layers
	std::string thg = "none";       ///< none | antenna | bulk : third-harmonic generation, not together with shg
	std::string update = "field";   ///< field | polarization : how the nonlinearity enters the update
	std::string chi_structure; ///< tensor names: a [chi2|chi3.<name>] section or a preset, config/chi_tensors.h; unset = zero
	std::string chi_first; ///< tensor names: a [chi2|chi3.<name>] section or a preset, config/chi_tensors.h; unset = zero
	std::string chi_middle; ///< tensor names: a [chi2|chi3.<name>] section or a preset, config/chi_tensors.h; unset = zero
	std::string chi_last; ///< tensor names: a [chi2|chi3.<name>] section or a preset, config/chi_tensors.h; unset = zero
	bool graphene_optimization = false;
};

/** [electrons]: the free-electron response. The keys of the hydrodynamic model are in [hydrodynamic]. */
struct ElectronsSection
{
	std::string model = "local";    ///< local | nonlocal | gnor | hydrodynamic
	double diffusion = 0.0;         ///< diffusion coefficient, nonlocal, gnor, hydrodynamic
	double fermi_velocity = 1.0;    ///< nonlocal, gnor, hydrodynamic
	bool operator==(const ElectronsSection&) const = default;
};

/// [hydrodynamic]: terms of the hydrodynamic model, [electrons] model = hydrodynamic
struct HydrodynamicSection
{
	bool pressure = true;           ///< quantum, Thomas-Fermi, pressure term
	bool convection = true;         ///< convective (v.grad)v term
	bool magnetic = true;           ///< magnetic Lorentz-force term
	bool split_nonlinear = false;   ///< separate linear and nonlinear DFTs, with multipole output
	bool operator==(const HydrodynamicSection&) const = default;
};

/** [temperature]: the electron temperature. two_temperature, pump and probe run the two-temperature
 *  model ([two_temperature]); pump and probe also read [pump_probe]; load reads a temperature map. */
struct TemperatureSection
{
	std::string model = "none";     ///< none | two_temperature | load | pump | probe
	bool operator==(const TemperatureSection&) const = default;
};

/// [two_temperature]: [temperature] model = two_temperature | pump | probe
struct TwoTemperatureSection
{
	bool save_frames = false;       ///< write electron temperature frames, model = two_temperature
	int rk_stages = 1; ///< Runge-Kutta stages
	bool operator==(const TwoTemperatureSection&) const = default;
};

/// [pump_probe]: [temperature] model = pump | probe
struct PumpProbeSection
{
	double delay = 0.0;             ///< [s]
	bool operator==(const PumpProbeSection&) const = default;
};

using ProbePoint = fdtd::monitors::ProbePoint; ///< one [probe.<name>] point, monitors/inputs.h

using MonitorBox = fdtd::monitors::MonitorBox; ///< one [monitor.<name>] section, monitors/inputs.h

/// one [dipole.<name>] section of a dipole source, [source] type = dipole
struct DipoleSpec
{
	std::string name;
	int x = 0; ///< [cells] inner domain, PML excluded
	int y = 0; ///< [cells] inner domain, PML excluded
	int z = 0; ///< [cells] inner domain, PML excluded
	double e_x = 0; ///< amplitude of the injected E components, times the pulse
	double e_y = 0; ///< amplitude of the injected E components, times the pulse
	double e_z = 1; ///< amplitude of the injected E components, times the pulse
};

using BlockSpec = fdtd::geometry::Block; ///< one [block.<name>] section, geometry/block.h

/// one [wave.<name>] section of a wave_superposition source
struct WaveSpec
{
	std::string name;
	double direction_x = 0; ///< unit propagation vector
	double direction_y = 1; ///< unit propagation vector
	double direction_z = 0; ///< unit propagation vector
	double e_x = 0; ///< E amplitude
	double e_y = 0; ///< E amplitude
	double e_z = 0; ///< E amplitude
	double h_x = 0; ///< H amplitude in vacuum units, scaled by n of the first medium
	double h_y = 0; ///< H amplitude in vacuum units, scaled by n of the first medium
	double h_z = 0; ///< H amplitude in vacuum units, scaled by n of the first medium
};

/** One key of a shape section: its name, the member of GeometrySection the builders read it from
 *  the shared pool the unnamed shapes still use, a unit and a description. `offset` is added to
 *  the value on the way to the pool, a rotation whose zero differs from the builder's. */
struct ShapeKey
{
	const char* key;
	double GeometrySection::*pool;
	const char* doc;
	double offset = 0;
};

/// A shape with its own section: [<name>] holds `keys`; docs/geometry/shapes/<name>.png shows them
struct ShapeSection
{
	const char* name;
	std::vector<ShapeKey> keys;
	void (*finish)(GeometrySection&) = nullptr; ///< for a key the builder reads through a formula, after the keys are in the pool
};

const std::vector<ShapeSection>& shape_sections();            ///< the shapes that have a section
const ShapeSection* find_shape_section(const std::string& shape); ///< nullptr: the shape reads [geometry]

struct SimConfig
{
	std::map<std::string, std::map<std::string, double>> shape_params; ///< [<shape>] key = value (shape_sections())
	RunSection run;
	PathsSection paths;
	GridSection grid;
	ParallelSection parallel;
	BoxesSection boxes;
	SourceSection source;
	PlaneWaveSection plane_wave;
	WaveSuperpositionSection wave_superposition;
	SignalSections signal;
	GaussianBeamSection gaussian_beam;
	ModeSection mode;
	SpectrumSection spectrum;
	OutputSection output;
	MultipolesSection multipoles;
	GeometrySection geometry;
	MaterialsSection materials;
	CustomMaterialSection custom_material;
	AnapoleSection anapole;
	NonlinearSection nonlinear;
	ElectronsSection electrons;
	HydrodynamicSection hydrodynamic;
	TemperatureSection temperature;
	TwoTemperatureSection two_temperature;
	PumpProbeSection pump_probe;
	std::vector<ProbePoint> probes; ///< one [probe.<name>] section each
	std::vector<MonitorBox> monitors; ///< one [monitor.<name>] section each
	std::vector<WaveSpec> waves; ///< one [wave.<name>] section each
	std::vector<BlockSpec> blocks; ///< one [block.<name>] section each, in file order
	std::vector<DipoleSpec> dipoles; ///< one [dipole.<name>] section each
	std::vector<ChiTensor> chi_tensors; ///< one [chi2.<name>] or [chi3.<name>] section each
};

SimConfig load_sim_config(const Config& ini);           ///< \throws syntax_error, validation_error
void validate_sim_config(const SimConfig& c);               ///< \throws validation_error
std::string to_ini(const SimConfig& c);                    ///< fully resolved configuration
std::string template_ini();                                ///< commented template with defaults

/** \brief Angular frequencies of the requested spectrum: fundamental band, then the harmonic band, if any. */
std::vector<double> spectrum_omegas(const SpectrumSection& s, double c0 = 299792458.0, double two_pi = 6.283185307179586);

int pulse_code(const std::string& name);           ///< legacy signal_type
/** \brief Start index of every rank along one axis, PML included.
 *  Equal split by default, the first `total % ranks` ranks get one extra cell; explicit
 *  `chunks` ("0 40 90") are used as given. No internal boundary may lie on a plane in
 *  `forbidden`, the TF/SF planes: an equal split moves such a boundary by one cell and
 *  describes the move in `note`; explicit chunks throw validation_error.
 *  Every rank needs at least 2 cells. */
std::vector<int> rank_starts(const char* axis, int total, int ranks, const std::string& chunks,
							const std::vector<int>& forbidden, std::string* note = nullptr);
/// TF/SF plane positions, global, PML included, along x, y, z; empty along periodic axes.
std::vector<int> tfsf_planes(const SimConfig& c, int axis);
/// [geometry] with the keys of the shape's own section written into the pool the builders read
GeometrySection resolved_geometry(const SimConfig& c);
int shape_code(const std::string& name);           ///< legacy antenna_type
int source_code(const std::string& name);          ///< legacy TFSF
int harmonic_mode(const std::string& name);        ///< none 0, antenna 1, bulk 2
int dft_plane_mask_bits(const OutputSection& o);       ///< legacy dft_plot bits: 1 yz, 2 xz, 4 xy, 8 magnetic
int harmonic_order(const std::string& name);       ///< none 1, second 2, third 3
int process_harmonic(const NonlinearSection& n);   ///< 2 with shg, 3 with thg, 1 otherwise

bool split_override(const std::string& spec, std::string& section, std::string& key, std::string& value);
/** \brief Parse INI text and layer the overrides, "section.key=value", and the directories on top.
 *  Precedence: file < overrides < input_dir / output_dir. \throws syntax_error, std::invalid_argument */
Config build_config(const std::string& ini_text, const std::string& name, const std::vector<std::string>& overrides,
						const std::string& input_dir = "", const std::string& output_dir = "");
std::vector<int> parse_int_list(const std::string& sect, const std::string& key, const std::string& raw);

} // namespace fdtd::config

#endif // __SIM_CONFIG_H__
