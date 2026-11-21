# Roadmap

State of the `driver-rework` branch and what is left. `[x]` = done, `[ ]` = to do.
Measurements are 1 rank, `MARCH= OPT=-O2`, unless noted. §1 and §2 describe the code and the
target; §3 onward is the record of the work, kept as it was written (names in old entries are the
names of that day).

## Status (2026-11-14)

The restructuring is at a break point. The solver is organised in modules with a layered
dependency order (`docs/dev/Architecture.md`), the core reads no application state, every function,
variable and configuration key follows one naming scheme, and the run is described by one INI file.

**Checks that every change passes** (`docs/dev/Tests.md`):
- bit-identical results against the previous build on the working suites (classic and Drude+2CP,
  sources, material models, periodic, monitors, custom boxes, near-to-far, SHG / THG, hydrodynamic
  split, two-temperature, multipoles, restart, the meshes of all shapes);
- the unit tests, `make regression` (7 cases), `make check-layers` (0 findings), and a build
  without warnings (`-Wall -Wextra`).
An intentional change of results is recorded where it is made, with its size; the list is in §7.

**In progress (2026-11-21):** `docs/dev/ReworkPlan.md`: the materials are definitions loaded into
slots and the reports have their new names and metadata, done; module ownership of the state is next.

**Open, in the order planned:**
1. Physics review of the hydrodynamic model: the x / y / z copies of its scans differ by 240 to 345
   lines each and can only be folded once the differences are classified (§5).
2. Matadata: the oscillator-fitting library as the source of the material-to-model mapping (§9).
3. Nonlinear models on dispersive materials (Varin's nonlinear Lorentz drive), and a review of the
   nonlinear processes with the material models (§9).
4. Smaller items: the remaining composite shapes ported to the geometry layers, other propagation
   axes, mode-source restart, Doxygen headers (§8, §9).

**Known defects** are in §7; the ones a user can meet: the two-temperature model with the
unsplit hydrodynamic model crashes, a magnetic Lorentz material as a layer crashes, and
restart is not available for the hydrodynamic, two-temperature, layered-wave and mode-source runs.

---

## 1. Where the code stands

| | |
|---|---|
| Source files (2026-11-14) | 83 `.cc`, 110 headers, about 29.7 k lines; the largest: `materials/hydro/hydrodynamics.cc` 2.1 k, `grid/decomposition.cc` 1.5 k, `config/sim_config.cc` 1.2 k, `materials/parameters.cc` 1.1 k, `app/setup.cc` 1.0 k |
| Layers | primitives (`core`, `io`), core library (`grid`, `fields`, `materials`, `boundary`, `update`, `sources`, `monitors`, `geometry`), application (`config`, `app`); `tools/check_layers.py` reports 0 upward dependencies |
| Shared state | no global state object. The application owns `app::RunConfig` (const after configure) and `app::RunState`; each core module receives a typed input from the application (`grid::Lattice`, `update::Spec`, `materials::Inputs` / `ModelSpec`, `geometry::Scene`, `sources::Spec`, `monitors::Spec`) and owns its own state (`inc/<module>/state.h`) |
| Functions | every function is declared in its module's header and namespace; sources and monitors are `Source` / `Monitor` objects in lists the driver builds |
| Naming | snake_case functions, variables and configuration keys; PascalCase types (`docs/dev/StyleGuide.md`) |

**What is still in the old style.** The hydrodynamic model (`materials/hydro`) keeps one copy of each
scan per field component and its own per-component array names; the composite shapes in
`geometry/composite.cc` and `legacy.cc` are point-wise builders with positional parameters, not yet
the set operations of `geometry/csg.h`; `app/setup.cc` and `materials/parameters.cc` are long
procedural set-up code. §2 is the target and §3 the order of work. Everything after §3 is history
and the remaining local items.

---

## 2. Target architecture (application layer first)

So far the design has been done bottom-up: kernels were cleaned without an agreed module structure above them. The target structure is fixed first, so each cleanup lands in its final place.

### 2.1 Modules

Each module is a header that documents the interface, plus an implementation file. A module owns its state; nothing is reached through `fdtd::state` once it is migrated.

| Module | Owns | Interface (sketch) | Replaces |
|---|---|---|---|
| `app` | CLI, `Driver`, `AppMachine`, signal state | phases `configure / decompose / init / step / finish`, `interrupted()` | `driver.cc`, `timeLoop.cc` wiring, `main.cc` |
| `config` | `SimConfig` (INI + overrides) | already a module | `pids` / `INI_DATA_STRUCT` (retired 2026-11-07; `applyConfig.cc` fills `app::RunConfig`) |
| `grid` | `Subdomain`, decomposition, halo exchange | `Subdomain::cpml(axis)` half-open ranges, `exchange(Field&)` | `decomposeDomain.cc`, `sendrecMPI.cc`, MPI datatypes in state |
| `fields` | E, H (and harmonic copies), history (`prev`, `prev_2`) | `Fields::E(axis)`, `rotate_history()` | `Ex`..`Hz`, `*_prev*`, `*_nl` globals |
| `materials` | material table, per-cell material maps, model registry | `MaterialMap`, `Registry::model(code)` | `materials.cc`, `matind_*`, per-model coefficient globals |
| `materials/models/*` | one header per model: point update now; coefficients and state to follow | `struct Drude { template <int A, class Boundary> static void update(const Cell &, Boundary &&); }` (done) | the per-model files (done 2026-11-07), `adeCoefficients.cc` (coefficients), the per-model state arrays |
| `nonlinear` | chi tables, nonlinear cell lists | `NonlinearTerm`, a P_NL/E-approach policy | `tensorNonlinear.cc`, `isotropicNonlinear.cc`, `chi*` globals |
| `boundary` | CPML profiles and psi state, periodic/Bloch | `Cpml::factors(axis)`, `Cpml::psi_terms(...)` | `cmpl.cc`, `boundaryTerms.cc`, CPML coefficient globals |
| `update` | the one E and H update path | `UpdateE(fields, materials, boundary, sources)` over stretches | `updateE.cc`, `updateH.cc`, `stretch_update.h` |
| `sources` | TF/SF plane wave (`TfsfSource`, `tfsf_h/e`, `AuxGrid`), superposition / Gaussian beam, layered-medium oblique, dipoles, mode source (all as `Source` objects, done); signals | `Source::inject_h/e(StepTime)` (done) | `planeWave*.cc`, `gaussianBeam.cc`, `setupSources.cc`, `signal.cc`, `OneDimensionalSimulator*` |
| `geometry` | layers, solids, shapes | already a module (`geometry.h`, `render_geometry.h`) | the remaining composites in `complexGeometries.cc` / `legacyGeometries.cc` / `readInGeometry.cc` |
| `monitors` | DFT planes/volumes, probes, flux boxes (Q, ART), multipoles, near-to-far | `Monitor::sample(step)`, `Monitor::finish()` | `updateDFT.cc`, `normalizeDFT.cc`, `qCoefficient.cc`, `art.cc`, `*Multipoles.cc`, `writeNearFar.cc`, `updateFieldMonitors.cc` |
| `io` | HDF5 writers, paths | already mostly a module (`h5io.h`) | `write*.cc` |
| `thermal` | two-temperature model | `TwoTemperature::step()` | `twoTemperature.cc`, `setupTemperature.cc` |

Layout (decided 2026-11-06): `inc/<module>/` for the headers, `src/<module>/` for the implementations. Everything lives in the project namespace `fdtd`; a module uses a sub-namespace (`fdtd::grid`, `fdtd::materials`, `fdtd::state`, ...) where it needs one. Calls are qualified by origin: `std::f()` for the C++ library including the wrapped C headers (`std::signal` from `<csignal>`), library prefixes as they are for MPI and HDF5 (`MPI_Allreduce`, `H5Fclose`), `::f()` only for a plain C function with no `std::` form, `fdtd::f()` (or the module's namespace) for ours. See `docs/dev/StyleGuide.md`.

### 2.2 One update path, models as policies

The update is written once: iterate stretches, compute the curl, call the material's point update, add the boundary terms. Everything that varies is a type or a callback:
- **Component:** `template <Axis A>` with a traits struct giving the arrays, the curl partners and the CPML factors. No X/Y/Z copies.
- **Material model:** a policy type with `Coeffs` (precomputed per material), `State` (its ADE arrays, ideally packed per stretch) and `point<A>(k, curl, boundary)`. Drude+2CP is already in this shape (`stretch::Drude2CPPoint` + `AdeRow`). Drude, Lorentz, modified Lorentz, Drude-Lorentz, magnetic Lorentz and the hydrodynamic model follow it.
- **Dispatch:** once per stretch, not per cell. The registry maps a material code to a model kind; a `switch` (or a `std::variant` visit) outside the k loop calls `UpdateStretch<Model, A>`. Each instantiation is a branch-free loop.
- **Boundary terms** (CPML psi, TF/SF corrections): callbacks given by `boundary` and `sources`, guarded per row and per k as now.
- **Linear/harmonic field sets:** the same update called with a different `Fields` view, not a copy.
- **Precision:** fields are `double` (the `REAL`/`REAL2` aliases are gone). Single precision or complex fields (Bloch, removed for now) come back as a template parameter of `Fields` and the update path, selected at compile time through the same mechanism as the models.

H, the DFT planes, the flux integrals and the probes follow the same idea: one loop, the variation as types or callbacks.

### 2.3 Graceful signals (`app`)

Done 2026-11-07 (`inc/app/signals.h`, `src/app/signals.cc`): tested by sending `SIGTERM` to the process on 1 rank (stopped at step 1705, exit 143) and to one of four ranks (all four stopped at step 4028, outputs written and readable, exit 143); `test_machine` covers the interrupted path.

- [x] **Handler:** for `SIGINT`/`SIGTERM` (and `SIGUSR1` for "write a checkpoint and continue", if wanted). Per [`std::signal`](https://en.cppreference.com/cpp/utility/program/signal) it only sets a `volatile std::sig_atomic_t` (or a lock-free `std::atomic<int>`); nothing else is safe there.
- [x] **State machine:** the step phase checks the flag at step boundaries. `AppMachine` gets an `interrupt` event and an `Interrupted` state that goes to `Finishing`. `Finishing` closes HDF5 files, writes what the monitors have (DFT normalised so far, probes), logs the step reached and exits with `128 + signal`.
- [x] **MPI:** a signal may reach only some ranks. Agree on it every step (one `MPI_Allreduce`; the loop already has two barriers per step) (`MPI_Allreduce` of the flag, MAX), so all ranks leave the loop at the same step and no rank blocks in a halo exchange. N small enough to react within a second or two.
- [x] **Test:** send `SIGTERM` mid-run on 1 and on 4 ranks. Check: no hang, outputs readable, exit code 143.

### 2.4 Transitional technique

Ownership moves module by module without touching every user at once:
- the module owns the data;
- `fdtd::state` keeps a reference alias (`inline auto &Ex = fields().E(X);`) until the last user is migrated;
- then the alias goes.

The `using namespace fdtd::state;` lines (one per translation unit, marked transitional) disappear in the same way: when a file uses no aliases, its directive goes.

---

### 2.5 Layers: the solver core as a library (user, 2026-11-13)

The solver is meant to be read as two layers, even though only one application is planned: the
**core** (an FDTD library: what a time step does) and the **application** (`fdtd simulate` and the
tools: what a run is). Other applications could be built on the core; designing for that keeps the
wiring explicit and shows where two pieces do the same job. Proposed layering, bottom to top:

| Layer | Modules | Knows about |
|---|---|---|
| 0. primitives | `core` (arrays, views, constants), `io/h5io`, `geometry/csg` | nothing else |
| 1. core library | `grid`, `fields`, `boundary`, `materials` (with its `nonlinear`, `hydro`, `thermal` parts), `update`, `sources`, `monitors`, `geometry` (builders) | layer 0 and each other, downward only (below) |
| 2. application | `config` (INI schema), `app` (driver, state machine, CLI, output tree, restart, logging policy), `tools/` | everything |

Rules for the core:
- **No upward includes.** A core file includes no `app/` or `config/` header. Today every core module
  reads `app::RunConfig` (the application's whole configuration, 257 fields, many of them derived per
  rank) and several read `app::RunState`; `fdtd::legacy` parameters are read in 35 places.
- **Parameters come in typed, per module.** Each module gets the part it needs as its own struct,
  built by the application from the configuration: e.g. `grid::Decomposition` (cells, PML, ranks),
  `boundary::CpmlSpec`, `sources::Signal` and `WaveFrame` (done), a `MaterialTable`
  (parameters and coefficients), `monitors::MonitorSpec` (planes, boxes, spectrum). `RunConfig` stays
  in the application as the parsed file and stops being the core's input.
- **Time and run state come as arguments** (`StepTime`, `SampleTime`, the step; done).
- **Side effects through the application:** file names, output directories and the log are given to
  the core (an output sink, a logger), not chosen inside it.
- **One order between core modules:** `grid` <- `fields` <- `materials` <- `boundary` <- `update`;
  `sources` and `monitors` on top of those; `geometry` writes a `MaterialMap`; the nonlinear,
  hydrodynamic and thermal parts belong to `materials`. A dependency against this order is a sign
  of a misplaced function.

What it should expose (redundancy candidates seen so far): `RunConfig` holds both input and derived
per-rank values (with `Subdomain`); the plot / ART / flux-box planes are computed in `grid` from
`geometry` and `analysis` settings, which belongs to `monitors`; the 1D auxiliary grid and the layered
1D simulators solve the same problem (§8); materials parameters are set by shape-dependent branches
in `materials/parameters.cc`, which is configuration, not core.

The check: `tools/check_layers.py` lists the include edges and the configuration reads per module
(a report now, a failing check once a module is clean). Done module by module, starting from the
bottom (`grid`, `fields`, `boundary`), each step bit-identical.

First report (`make check-layers`, 2026-11-13, after `globals.h` was removed):

| module | app includes | app reads | against the order |
|---|---|---|---|
| io | 6 | 1 | 0 |
| grid | 4 | 10 | 4 |
| fields | 0 | 0 | 0 |
| materials | 5 | 34 | 14 |
| boundary | 1 | 2 | 3 |
| update | 1 | 3 | 4 |
| sources | 17 | 30 | 0 |
| monitors | 19 | 29 | 0 |
| geometry | 7 | 6 | 0 |
| nonlinear / hydro / thermal | 2 / 3 / 5 | 3 / 4 / 6 | 0 |

What it showed beyond the `RunConfig` reads, and what was done (2026-11-13, second batch):
- [x] `io/output_layout` (where the rank's cells sit in the output files) is a grid fact: now
  `grid::output_layout()` (`inc/grid/output_layout.h`, `src/grid/lattice.cc`).
- [x] The model kernels included `update/stretch.h`: the row-loop helper is the models' own, now
  `inc/materials/stretch.h` (namespace `fdtd::materials::stretch`). `update` sits above `materials`.
- [x] `fdtd::legacy` removed (`app/legacy_parameters.h`): the setup failure codes duplicated the error
  message every failure already gives (`Agree(ok, why)`), so they are gone; `grid::set_lattice_data`
  returns the reason as text. The per-rank cell tables are vectors (no 512-rank limit); the one real
  capacity limit, the per-frequency tables, is `fdtd::max_frequencies` (`inc/core/limits.h`). The
  unused 131072-rank limit is gone.
- [x] `grid::Lattice` (`inc/grid/lattice.h`): cells, CPML, TF/SF inset, periodicity, cell size and its
  inverse, set by the driver once after `setup_domain` and read by `grid` (halo exchange, output
  layout, decomposition), `boundary` and `update` instead of `run_cfg.grid`. The first module input
  in the §2.5 sense. `RunConfig::grid.num_cells_total_*` (a copy of `nr_cells_tot_*` made after the
  decomposition) merged into `nr_cells_tot_*`.
- [x] The decomposition read the geometry and analysis settings (91 reads) to place the planes: the
  placement (TF/SF box, ART planes, flux boxes, plot planes, some following the antenna shape) is
  application policy, now `app::global_planes()` (`src/app/planes.cc`, same arithmetic); the grid gets
  the result as `grid::GlobalPlanes` (`inc/grid/planes.h`) and only maps it onto the ranks.
  `grid` reads no configuration.
- [x] (eighth batch) `io` (a primitive) still includes `materials/material_field.h` (the archive serialises material
  fields) and `grid/state.h` (the rank number in log lines): both should be passed in.
- [x] `grid/halo_exchange.cc` included `fields`: the exchanges of the fields object are
  `fields::exchange_halos` / `exchange_longitudinal_e` (`src/fields/halos.cc`) over the grid's
  `exchange` primitive. `grid` and `fields` are clean; `make check-layers` fails if they regress
  (`--strict core grid fields`).
- [x] (ninth / tenth batch) `boundary` and `update` include `sources/aux_grid.h` and `sources/tfsf.h`: the TF/SF corrections
  are called from inside the boundary and update loops instead of being a separate step of the time loop.
- [x] `update` and `boundary` read physics switches: `update::Spec` (`inc/update/update.h`: TF/SF,
  SHG/THG, P_NL mode, magnetic Lorentz, harmonic fields, cellwise materials, background media,
  `coeff_h`) and `Cpml::background` / `Cpml::harmonic`, set by the driver after the run is
  configured. Neither module reads the configuration; what is left are the include edges to
  `sources` (TF/SF, the 1D grid's layer materials) and `nonlinear`.
- [x] (sixth / ninth batch) `materials/component.h` includes `app/run_config.h`, `boundary/cpml.h` and `hydro/state.h`: the
  per-component view gathers everything a kernel touches; it is the natural place for the typed
  parameters (one struct per model) that replace those reads.

Report after the second batch (`make check-layers`; field reads = `run_cfg.<group>.<field>`):

| module | app includes | field reads | against the order |
|---|---|---|---|
| io | 2 | 0 | 0 |
| grid | 1 | 91 | 2 |
| boundary | 0 | 4 | 2 |
| update | 0 | 25 | 4 |
| materials | 4 | 320 | 10 |
| sources | 16 | 269 | 0 |
| monitors | 18 | 369 | 0 |
| geometry | 6 | 1394 | 0 |
| hydro / thermal / nonlinear | 2 / 4 / 1 | 281 / 31 / 8 | 0 |
| total | 54 | 2792 | 18 |

Report after the third batch:

| module | app includes | field reads | against the order |
|---|---|---|---|
| io | 2 | 0 | 0 |
| grid, fields | 0 | 0 | 0 |
| boundary | 0 | 0 | 2 |
| update | 0 | 0 | 4 |
| materials | 4 | 320 | 10 |
| sources | 16 | 269 | 0 |
| monitors | 18 | 369 | 0 |
| geometry | 6 | 1394 | 0 |
| hydro / thermal / nonlinear | 2 / 4 / 1 | 281 / 31 / 8 | 0 |
| total | 53 | 2672 | 16 |

After the fourth batch (2026-11-13): the removed shapes took ~500 configuration reads with them, and
every core module reads the grid facts (cell size and inverse, CPML, cells, TF/SF inset, periodicity,
inner cells `Lattice::inner`) from `grid::lattice()` instead of `run_cfg.grid` (24 files). Field reads
1513 (from 2672): geometry 562, monitors 325, materials 306, sources 168, hydro 113, thermal 31,
nonlinear 8.

Fifth batch (2026-11-13): `geometry` reads no configuration and is strict in `make check-layers`.
- The shape dispatch (`create_structures`: shape and dimensions -> builder calls, the nonlinear
  antenna-material rule) is application policy: `src/app/structures.cc` (was `geometry/read_in.cc`).
- `geometry::Scene` (`inc/geometry/scene.h`): shape, draw flag, staircasing, background media,
  centre cell and the round_hybrid gap material, set by the application before drawing; the builders
  read it instead of `run_cfg`.
- The `[anapole]` parameters and the `[block.*]` boxes are geometry types (`geom::AnapoleSpec`,
  `geometry::Block`) that the configuration schema uses (`using AnapoleSection = ...`); the builders
  take them as arguments and the application holds them (`app::anapole_input`, `app::blocks`).
- Removed as unused: `read_geometries` (`Geometries.txt`), `read_2d_geometry` (`2D_Geometry.txt`),
  `create_fcc_sphere_cluster`.
- Field reads 951: monitors 325, materials 306, sources 168, hydro 113, thermal 31, nonlinear 8.
- Seen on the way: `nano_antenna` is `shape != none`, and the antenna builders return without
  drawing when it is false; `[block.*]` boxes are drawn either way (they do not check it).

Sixth batch (2026-11-13): `materials` reads no configuration.
- `materials::Inputs` (`inc/materials/inputs.h`): the media, the dielectric indices, the user values
  of the built-in models (Fermi velocity, GNOR diffusion, ITO carrier densities, graphene stretch,
  the middle-layer thickness that selects the ITO film model), `[custom_material]` and the chi2 / chi3
  tables of the four media. `set_parameters(const Inputs&)`; the application builds it.
- `materials::ModelSpec` (`inc/materials/spec.h`): time step, Fermi velocity, SHG/THG, harmonic
  fields, cellwise materials, nonlocal / hydrodynamic / GNOR, two-temperature / loaded temperature;
  the model kernels and the coefficients read it. Set after the time step in `initialise_run`.
- `read_magnetisation_field(input_dir)` takes the folder.
- The TF/SF H correction coefficients (`chx/chy/chz`, used only by the TF/SF source) were set by the
  material coefficients into `update`: now `sources::ch*`, set by the application.
- Two switches that were always false removed with their branches: `e_dens_ft` (set from the
  configuration, then forced false in `initialise_run` before anything read it) and
  `no_interband_transitions` (forced false the same way). With them the only call of
  `hydro::write_lwa_multipoles` went (it needed `e_dens_ft`): the LWA multipoles were unreachable,
  which belongs to the hydro review (§3).
- [x] The nonlinear, hydrodynamic and thermal parts are part of `materials` (user, 2026-11-13: they
  enable analysis of material properties and sit beside the choice of material model, Drude, Lorentz
  and the others): `inc/materials/{nonlinear,hydro,thermal}/`, `src/materials/{...}/`, namespaces
  `fdtd::materials::nonlinear`, `::hydro`, `::thermal`. The order findings between them and the
  material kernels are gone; what `materials` still reads (136 field reads) is in these parts
  (`hydro/hydrodynamics.cc` 95, `thermal/two_temperature.cc` 17, ...), the per-component copies
  waiting for the physics review (§5).

Seventh batch (2026-11-13): `sources` reads no configuration and no run state; strict in
`make check-layers`.
- `sources::Spec` (`inc/sources/spec.h`): time step, H coefficient, 1D grid steps and cells, the
  three media, middle layer, centre, the band (SHG/THG, frequency count, fundamental, highest,
  CW turn-on), which source is on and its parameters (TF/SF, mode, dipoles, Gaussian beam,
  superposition, oblique through the layers, circular delay, the CW flags, the beam samples), the
  plot periods, the report and input folders and the log start. Set at the end of `initialise_run`.
- The signal: `Signal(Params)` builds it (shape, dt, peak, CW frequency and turn-on, chirp, CW
  length, band edge, Gaussian width and carrier); `app::configure_signal()` (`src/app/signal_config.cc`)
  fills the parameters from the configuration and sets the run's signal. The CW envelope read the
  fundamental and the ramp read the peak from the configuration although they were given them:
  they use their arguments now (the same values).
- `transform_pulse(steps)` returns the step count (a CW run sets its own length) instead of writing
  the run state; `build_sources(..., steps)` gives the dipole signal table its length.
- `[mode]` is a sources type (`sources::ModeSpec`, `inc/sources/mode_spec.h`) that the schema uses.
- Field reads in the core: 455, of which monitors 319 and the materials parts 136.

Eighth batch (2026-11-13): `monitors` and `io` read no configuration and no run state; `make
check-layers` is strict for core, io, grid, fields, geometry, sources and monitors.
- `monitors::Spec` (`inc/monitors/spec.h`): time step and step count, the band, the analyses asked
  for (absorption, ART, sphere multipoles, DFT planes, field output), the flux-box inset, the shape
  and structure (centre, shift, thickness, materials), the background media, the source peak,
  frequency and polarisation, the output block size and the field / frame / report folders. Set in
  `Driver::configure_run` once the pulse transform has fixed the step count.
- `[probe.*]` and `[monitor.*]` are monitors types (`monitors/inputs.h`) that the schema uses.
- The file-name scratch buffer of the run state (`run_state().cstring`) became a buffer of the two
  monitors that format step numbers.
- `io`: the log reads the rank from MPI instead of `grid::world_rank`; the restart archive
  (`app/archive.h`, was `io/archive.h`) is the application's.
- Redundancy seen: `RunConfig::grid.nr_cells_no_cmpl_*`, `inner_dim[]` are the lattice's inner cells
  (`Lattice::inner`); `inner_dim_nf[]` is that minus the flux-box inset; `globalsizes[]` is inner + 1
  on ranks that are not CPML-only. The monitors now use the lattice for the first two; the copies
  remain for the application and the thermal output.
- What is left in the core: 136 field reads in the hydrodynamic, thermal and nonlinear parts of
  `materials`, and 8 includes against the order (component.h and hydrodynamics.cc include the CPML;
  boundary and update include the TF/SF and 1D-grid headers of `sources`).

Ninth batch (2026-11-13): no core module reads the run configuration or the run state (0 of the
2792 field reads of the first report); 2 includes against the order remain. Strict in
`make check-layers`: core, io, grid, fields, update, geometry, sources, monitors.
- The hydrodynamic, thermal and nonlinear parts read `materials::ModelSpec` (extended with the
  pump-probe flags, the structure and background material codes, staircasing and the TF/SF y
  shrink). The thermal output computes its block size from the lattice as the application did
  (inner + 1, 0 on CPML-only ranks).
- The background along y (`mat_y_uni`, held by the 1D plane-wave grid) is a materials array,
  `materials::background_y`; the 1D grid keeps a reference to it. Boundary, update, geometry and
  hydro read it from materials.
- The 1D grid steps from the time loop (`Driver::update_magnetic` / `update_electric`, first, as it
  did at the start of `update_h` / `update_e`); the TF/SF polarisation flags are in `update::Spec`.
  `update` includes nothing of `sources`.
- Left, both real couplings to decide:
  1. `boundary/ade_terms.cc` adds the TF/SF E correction (`sources::tfsf_e_terms`) inside the CPML
     row term, so the boundary reads the incident field. The TF/SF surface is a boundary of the
     scattered-field region: one way is a `boundary::IncidentField` interface the TF/SF source
     implements and the application hands to the boundary, or the TF/SF correction moves into
     `boundary` with only the 1D grid and the signal staying in `sources`.
  2. `materials/hydro/hydrodynamics.cc` uses the CPML stretching (`inv_kappa`) in its derivatives.
     The stretching profile is a property of the coordinates, not of the absorbing fields: the
     profile (kappa, b, c) could live below materials (in `grid`) and the CPML keep the psi fields.

Tenth batch (2026-11-13, both decided as proposed): `make check-layers` finds nothing and is strict
for every core module.
- The TF/SF correction reaches the boundary as its incident field: `boundary::IncidentETerms`
  (`inc/boundary/incident_field.h`), a function the application hands over
  (`set_incident_e_terms(&sources::tfsf_e_terms)`); `cpml_tfsf_terms` adds it after the CPML psi.
  The boundary no longer includes `sources`.
- The coordinate stretching (1 / kappa per axis, face and grading position) is `grid::Stretching`
  (`inc/grid/stretching.h`), filled by `Cpml::set_profiles`; the CPML profile keeps b and c (the
  psi). The field updates, the hydrodynamic model and the 1D grid read it from `grid`.
- Removed: the allocation of every CPML profile table in `allocate_fields`, which
  `set_profiles` replaced right after.

`RunConfig` after the layering (2026-11-13): fields that no code read are gone, as are the copies of
lattice facts: `grid.nr_cells_no_cmpl_*`, `inner_dim[]`, `inner_dim_nf[]`, `globalsizes[]`
(`Lattice::inner`; the output block size is inner + 1 on ranks that are not CPML-only, computed
where it is used), `grid.inv_d*` (`Lattice::inv_d`), `source.raised_cosine` (written, never read),
`source.PERIOD_BURST`, `analysis.DFT`. 164 fields remain, read by the application only.
- Field reads 643: monitors 325, sources 168, hydro 111, thermal 31, nonlinear 8.

## 3. Order of work

1. [x] **Module layout** (2026-11-07): the Makefile builds `src/*/*.cc`. First module `boundary` (`inc/boundary/cpml.h`, `src/boundary/{cpml,e_faces,ade_terms}.cc`, namespace `fdtd::boundary`): `Cpml` owns the grading profiles of every axis and face and all psi arrays (linear and harmonic); one grading routine replaces 8 copied loops, the allocation is generated from the face/field/component pattern (48 hand-written lines), `release()` replaces 46 lines. The 70 legacy names are `inline` references into `cpml()` in `fdtd::state` (the transitional technique of §2.4). Bit-identical.
2. [x] **Signals** (§2.3): `fdtd::app`, the `Interrupted` state.
3. [~] **`fields` and `grid`** (2026-11-07):
   - [x] `fields` module (`inc/fields/fields.h`): `Fields` owns E, H (linear and harmonic sets), E_prev, E_prev_2, H_prev and the E scratch copies; `allocate(needs)` replaces 27 scattered allocations. The 27 legacy names are references into `fields()`.
   - [x] `grid` module: `inc/grid/subdomain.h`, `src/grid/{subdomain,decomposition,halo_exchange}.cc`; `grid::Subdomain`, `grid::this_rank()`, `grid::member_comm`.
   - [x] CPML `min`/`max` are the true local interface indices on every rank that has inner cells along the axis; ranks that are all CPML along it keep an empty range (with the true index, their interface cell passed the inner test and the E scan read a material map they do not have: segfault, found by the test below). Checked with rank edges exactly on the interfaces (`test/edge_layouts.sh`: x, y, z and all three, vacuum and Drude background, SHG, hydrodynamic split) and the usual layouts, all bit-identical to 1 rank.
   - [x] Halo exchange: the six `SendRecv{E,H}on{X,Y,Z}` copies (330 lines) are one `grid::exchange(Field, axis, ...)` (`inc/grid/halo_exchange.h`). Checked against the previous build with periodic x and z on 111, 211, 311, 112, 113, 212.
   - [ ] `MultiSwitch` / `MultiSwitch2` (the margin-grid exchange of the hydrodynamic model and the LWA multipoles, about 430 lines): same idea, but the branches are irregular (two or three planes, different order per component); fold it with the hydrodynamic review.
   - [x] Halo exchange on the fields object (2026-11-13): `grid::exchange_halos(Field, set)` (the two transverse components on each axis, for the linear or harmonic set) and `grid::exchange_longitudinal_e()` (the nonlinear layout); the time loop calls these instead of fifteen component-wise calls.
   - [x] The `_extra` values named (2026-11-13): `*_cpml_min_extra` / `*_cpml_max_extra` are `*_nea_cpml_end` (last local index of the near CPML) and `*_far_cpml_start` (first index of the far CPML, also the origin of the far psi arrays).
   - [x] `set_lattice_data` split (2026-11-13): the 1100-line function is `LatticeSetup` with eight steps (`rank_extent`, `global_planes`, `near_to_far`, `tfsf_box`, `local_planes`, `plane_wave_faces`, `scattering_box`, `total_box`) sharing the rank's extents; bit-identical on the regression, source, monitor, near-to-far and periodic suites. The per-axis copies inside the steps remain.
   - [ ] Half-open ranges at every call site (the E scans use `[min, max]`, the nonlinear loops `[min, max)`): left as is; changing a bound is a numerical change per loop, to do with the module that owns the loop.
4. [~] **`materials` + `update`:**
   - [x] One E scan (2026-11-07): `ex_update_scan`/`ey_update_scan`/`ez_update_scan` (three ~130-line copies) are `e_scan<A>` with an `EAxis<A>` traits struct (arrays, model entry points, boundary-term argument order). The rule that made them differ: along A the loop is `0..cells` with an exclusive inner test, across A `1..owned` with an inclusive one; the curl is `dH_{A+2}/d(A+1) - dH_{A+1}/d(A+2)` with `Diff<axis>` resolved at compile time. Bit-identical for every model as structure and background, nonlocal/GNOR, Drude background at 0/45/90 deg and periodic, SHG/THG matrix, dielectric exact SHG/THG, hydro split, edge layouts.
   - Performance pitfall found: the merged template first ran 17% slower because GCC stopped inlining the per-cell kernels into the larger body (gprof: 155 M calls each to the derivative lambda and to `Drude2CPPoint`). The hot curl now calls `Diff<axis>` directly and the stretch kernels are `__attribute__((always_inline))`; timing is back to the previous build (vacuum and Drude background).
   - [x] Files moved (2026-11-07): `src/update/update_e.cc`, `src/update/update_h.cc`, `inc/update/stretch.h`.
   - [x] Compiler hints written directly (`__attribute__((always_inline)) inline`, `[[assume]]` after an `assert`); style guide "Compiler hints".
   - [x] Models as policy types (2026-11-07): `inc/materials/models/{drude,drude_2cp,lorentz,modified_lorentz,hydrodynamic,drude_lorentz,magnetic_lorentz}.h`, each `struct Model { template <int A, class Boundary> static void update(const Cell &, Boundary &&); }` (shape in `models/cell.h`). The arrays come from `materials::Component<A>` (`inc/materials/component.h`: one generic struct, `pick<A>(x, y, z)`, no per-axis specialisations); the per-component parts are in cyclic form: the magnetisation cross term is `M_{B1} dE_{B2} - M_{B2} dE_{B1}`, the nonlocal grad(div P) one stencil over the axes. `materials::with_model(mm, f)` (`inc/materials/models.h`) maps a material index to its model; `e_scan` calls it once per stretch and loops the cells inside, so each model's loop is its own instantiation. The 21 `Update*_X/Y/Z` functions, the six model `.cc` files and the model lines of `EAxis` are gone. Bit-identical: every model as structure (100-400) and background (100-350), nonlocal, GNOR, hydro split, SHG/THG matrix, regression cases; a Lorentz background runs at the same speed.
   - Boundary terms per model: Drude+2CP, modified Lorentz, hydrodynamic and nonlinear Drude-Lorentz add the CPML psi and TF/SF terms inside their update (the boundary callback); Drude and Lorentz get them after the scan, from the E-face psi pass (`boundary::e_face`, scaled by `ce_2_d`) and the TF/SF faces (`sources::tfsf_e`); magnetic Lorentz gets neither (§7). One route for all models comes with the boundary callback in the stretch loop.
   - [x] Coefficients into the models (2026-11-07): every model has `struct Coeffs` and `static Coeffs compute(mm, code, dt)` in its header, with the expressions that were in `compute_ade_coefficients`; `materials::coeff_table<Model>` (`inc/materials/coefficients.h`) holds one `Coeffs` per material index and `coeffs<Model>(mm)` reads it. Shared shapes: `Poles` (b0..b4 per pole, c0..c3) for Drude+2CP, hydrodynamic, Drude-Lorentz and magnetic Lorentz; `ElectronGas` (nl1..nl5) for the hydrodynamic and Drude-Lorentz models (`materials::electron_gas(mm)` for the hydrodynamics code); `DispersiveE` (ce1, ce2) for Drude, Lorentz and magnetic Lorentz. 33 global coefficient arrays are gone (`ALPHA`, `BETA`, `ce_1_d`, `ce_2_d`, `a_cp_*`, `b_cp_*`, `c_cp_*`, `aj`, `bj`, `alpha_j_*`, `beta_j`, `ce_*_l`, `nonlin_drude_*`); `compute_ade_coefficients` is a dispatch over the code ranges (268 -> 60 lines). The E-face pass, TF/SF faces, nonlinear term, 1D auxiliary grid and the Drude+2CP stretch kernel (which now takes the coefficients as an argument) read the model tables. Bit-identical (before the bug fix below): material models as structure and background, polarisation in dispersive backgrounds, SHG/THG matrix, regression cases, hydro split; 7% faster on a Drude+2CP background and sphere (one coefficient struct per material instead of a dozen arrays).
   - [x] Bug found on the way and fixed (2026-11-07, §7): the Drude+2CP pole sums were shared across materials.
   - [x] Dispersive state limited to the material (2026-11-07): the polarisations and critical-point poles of the point-update models (Drude+2CP, Drude, Lorentz, modified Lorentz, without the nonlocal / hydrodynamic / two-temperature arrays) are stored only for the rows (i, j) of the material's bounding box on the rank, full along k. `carray::boxed` (rank 3 and 4) keeps the indexing of a full array: the rows outside the box point to one shared row of zeros, which is what the full arrays held there (only a material's own cells are written). `material_field` stores one array per material (boxed or full). The allocation moved after the geometry (`app::allocate_material_state`, from `Driver::build_geometry`), which computes the boxes from the inner material maps; a material that can be a background (first / middle / second medium, the 1D layer profile) keeps full rows. `app::material_state_intact()` checks at release that no write landed outside a box (none in any test run). Rank 0 logs the memory: a 5 nm sphere in 36^3 0.8 MB instead of 20.8 MB, a 30 nm sphere in 100^3 58.8 MB instead of 221.8 MB; a Drude+2CP background unchanged. Bit-identical: regression cases, hydro split, material models, polarisation, SHG/THG matrix, periodic suite, all sources, monitors, custom boxes (two-sweep DFG), CPML-edge layouts with a 20 nm sphere vs one rank, unit tests (`test_carray` covers `boxed`). The sphere run is about 20 % faster (8.2 -> 6.4 s, 100^3, 1 rank). Full per-stretch packing was not done: the hydrodynamic and nonlocal code reads P at neighbouring cells, and the box gives most of the memory.
   - [x] The 1D auxiliary grid uses the models (2026-11-07): Drude, Lorentz and modified Lorentz have a `point(cf, ..., k, AdeRow, curl[, boundary])` kernel on a row and an index, like `stretch::Drude2CPPoint`; the 3D updates and `sources::AuxGrid` (the 1D line is one row indexed by j) both call it, so the four model copies in `aux_grid.cc` are gone. Bit-identical: regression cases, material models (Drude, Lorentz, modified Lorentz and Drude+2CP backgrounds through the 1D grid), polarisation, layered media, all sources; timing equal (Drude background).
   - [x] Row loops for Drude and Lorentz (2026-11-07): a model may have `stretch<A>(mm, i, j, k0, k1, curl)`; the stretch dispatch in `e_scan` calls it when it exists (`if constexpr (requires ...)`), else the per-cell update. Drude and Lorentz (no boundary callback) set their row up once per stretch and call `point()` per cell, like the Drude+2CP path. Bit-identical (regression cases, hydro split, material models, polarisation); Drude and Lorentz backgrounds 10-15 % faster (100^3, 1 rank: 10.3 -> 9.0 s and 10.3 -> 9.2 s). The dispatch of §2.2 (model chosen once per stretch, a branch-free loop per model) now holds for every model with a row loop; modified Lorentz, hydrodynamic, Drude-Lorentz and magnetic Lorentz keep the per-cell loop (their boundary callback, or stencils across cells).
   - [x] Kernels in the model types (2026-11-13): `models::Dielectric` (`applies`, `stretch`) and `models::Drude2CP` (`stretch_index`, `point`, `stretch`) hold what was in `inc/update/stretch.h`, which keeps only the stretch infrastructure (`for_each_stretch`, `StretchList`, `AdeRow`, `Psi`). Modified Lorentz runs on the row path with the boundary callback guarded per row and per k, like Drude+2CP; the E scan dispatches every model the same way (row loop with or without boundary callback, else per cell). Bit-identical (regression, every material model as structure and background, hydrodynamic split). Speed of the modified Lorentz row path not measured.
   - [x] Material maps (2026-11-13): `materials::MaterialMap` (`inc/materials/material_map.h`) owns the inner, cell-centre and full maps with `allocate_inner`, `allocate_full`, `extend_into_pml` (was in the mode source), `dispersive_boxes` (was in the state allocation) and `release`; the solver names (`matind_x` ...) are references into it.
   - [x] The dispersive state (P, J, poles, parameters) as one object: see `materials::DispersiveState` below.
5. [~] **`sources`, then `monitors`:** both are lists of objects with `inject`/`sample` hooks, which is already how `timeLoop.cc` assembles them as lambdas.
   - [x] The 1D auxiliary grid of the TF/SF source (2026-11-07): `inc/sources/aux_grid.h`, `src/sources/aux_grid.cc`, `fdtd::sources::AuxGrid` owns the two lines (source line, periodic reference line) and the 1D material map; `AuxLine` holds both polarisations as `[2]` arrays. One E and one H cell update serve both polarisations (they differ only in sign, `add<P>`) and both lines; before, four copies in `ElectricFieldUpdateScan` (about 520 lines) and two in `MagneticFieldUpdateScan` (about 90). `allocate`/`release` replace about 250 lines in `setupSimulation.cc` and `memoryRelease.cc`. The 45 legacy names are references into `aux_grid()`. Bit-identical: every model as structure and background, periodic x+z in vacuum and in Drude+2CP, Drude, Lorentz and modified Lorentz backgrounds, layered first/last media, regression cases, SHG/THG matrix, hydro split, layouts 222/411/141/114.
   - [x] TF/SF face corrections (2026-11-07): `inc/sources/tfsf.h`, `src/sources/tfsf.cc`, `sources::tfsf_h(Face)` / `tfsf_e(Face)`. One plane loop for all faces; the E coefficient of the background is one rule (`ce_2`, `ce_2_d` for Drude/Lorentz, none for models that add the terms themselves), the same as the E-face CPML pass. Replaces the twelve `*PlaneWaveExcitation{X,Y,Z}{nea,far}` functions (`planeWave.cc`, 697 lines, removed) and their dead Bloch branches. Bit-identical: polarisation 0/45/90 in vacuum, Drude, Lorentz, Drude+2CP and layered backgrounds, periodic x+z, every model, regression cases, SHG/THG matrix, hydro split, layouts 222/311/131/113.
   - [x] `sources::TfsfSource` (2026-11-07, `inc/sources/tfsf_source.h`): built once per sweep for this rank; owns the face lists and the injection of the signal into the 1D grid (both lines, oblique/circular polarisation); `inject_h(t)` / `inject_e()`. `BuildInjectors` pushes two calls instead of about 60 lines of lambdas.
   - [x] Plane-wave superposition and Gaussian beam (2026-11-07): `inc/sources/superposition.h`, `src/sources/superposition.cc`, `sources::superposition_h/e(t, waves)`. One H and one E face update for the six faces (normal n, tangential a < b, a sign per face); the delays are summed in x, y, z order as before. Replaces the twelve face functions of `gaussianBeam.cc` (819 lines, removed); the wave list is no longer copied every step. Bit-identical: three-wave superposition in vacuum and Drude background, with the derivative pulse and dispersion correction, layouts 222/131/113; Gaussian beam (with focal and waist factors set, §7).
   - [x] Superposition through a layered background (2026-11-07): `inc/sources/oblique.h`, `src/sources/oblique.cc`, `sources::oblique_h/e(ods, t)`, the incident field from `OneDimensionalSimulator`. The face table, signs and box geometry are shared with the superposition in `src/sources/faces.h` (module-internal). Replaces the twelve functions of `obliquePlaneWave.cc` (417 lines, removed). Bit-identical: one and three waves over a glass layer (last medium, and middle + last), layouts 211/121/112/222.
   - [x] One source interface (2026-11-07): `inc/sources/source.h`, `sources::Source` with `inject_h/inject_e(const StepTime &)`, and `build_sources(ods)` builds this rank's list once per sweep: `TfsfSource`, the mode source, the superposition / Gaussian beam, the layered-medium waves, the dipoles. The time loop calls the hooks on the list; the `std::function` lists and the source lambdas in `timeLoop.cc` are gone. Bit-identical: two dipoles (layouts 211, 121), superposition, layered waves, oblique polarisation, periodic Drude background, regression cases.
   - [x] The CPML + TF/SF terms of the ADE models (2026-11-07): `boundary::cpml_tfsf_terms<A>(i, j, k, c, p, q, e)` replaces the three per-component copies; the CPML psi of one derivative axis is one function (`Psi<A, D, sign>`: + for the B1 faces, - for B2, summed in the previous order), the TF/SF terms are `sources::tfsf_e_terms`. `EAxis` loses its `terms` lines. Pitfall: `cpml()` is out of line and was called per CPML cell (4% slower on a Drude+2CP background); cached once per file, timing equal. Bit-identical: every model, Drude+2CP / modified Lorentz / Drude-Lorentz backgrounds at 0/45/90 deg, periodic, edge layouts in Drude+2CP and modified Lorentz, SHG/THG matrix, hydro split.
   - [x] Dipoles own their data (2026-11-07): positions and amplitudes from `dipole_input`, the signal tabulated per step when the source is built; `SetupDipoleSources`, the five globals and the `[step][dipole][6]` table are gone.
   - [x] The wave list (2026-11-07): `sources::build_waves` (`inc/sources/waves.h`) builds the superposition / Gaussian-beam waves, the focus, the TF/SF box positions and the beam offset, moved out of `Inizialization()` unchanged; the list lives in the driver context next to the 1D simulations and is passed to `build_sources`. The `PLANE_WAVES` global and nine file-local TF/SF indices are gone.
   - [x] Sources (2026-11-13): the signal is a `Signal` value (`inc/sources/signal.h`: shape, timing, `value(t)`, `derivative(t)`, `extent`, `fits`), built once by `configure_signal()` after the time step is set; every source holds it (`build_sources(signal, waves, frame, ods)`), and the 15 signal globals are gone (`run_signal()` serves the monitors' normalisation and the tools). The superposition anchor (focus, TF/SF faces, common delay) is a `WaveFrame` held by the superposition source and by the 1D simulators, filled by `build_waves(waves, frame)`; the 10 frame globals are gone.
   - [x] `monitors`, first part (2026-11-07): the flux boxes (`inc/monitors/flux_boxes.h`, `src/monitors/flux_boxes.cc`). `monitors::FluxBoxes` owns the DFT planes of the scattering and total boxes: `box[b].face[normal][side].re/im[E/H][field/harmonic][component]`. One `allocate`, `release` and `sample` replace about 190 lines of allocation, 125 of release (which missed the harmonic planes) and the 48 written-out plane-DFT calls of `DFT_Update`. The 192 legacy names are references into it. Bit-identical: the Q coefficient files of every regression case and of the hydrodynamic split (both sets), layouts 222/311/113, all sources.
   - [x] The Q integrals (2026-11-07): one `face_flux(set, box, normal, side)` over the `FluxBoxes` faces, `Re(E_a H_b* - E_b H_a*)` with the cyclic axes a = n+1, b = n+2, replaces `QcoefficientYZ/XZ/XY` and their 30 calls in `Q_coefficients` and `ART_full_TFSF`; the 52 per-face result arrays are one `q_face[set][box][normal][side]`. `qCoefficient.cc` 703 -> 413 lines. Every `Qcoeff*`, `ART-TFSF*` and `NL_Efficiency*` file byte-identical (regression cases, hydrodynamic split with its harmonic files).
   - [x] Display DFT planes (2026-11-07): `monitors::DftPlanes` (`inc/monitors/dft_planes.h`) owns the three animation planes (`plane[normal].re/im[E/H][field/harmonic][component]`) and the E-only far-field xz planes and graphene plane. One `allocate` / `release` / `sample` replace 90 written-out allocations, the release lines (36 of the 90 were never released) and 15 plane-DFT calls; the allocation function keeps only the animation flags (its P_NL allocation was a duplicate of the one with the volume DFT, on every rank, and is removed). The 90 legacy names are references into it. Bit-identical: all three planes with E and H, the harmonic planes of the hydrodynamic split, the forward/backward planes (cylinder), regression cases, all sources, layouts 222/141.
   - [x] ART planes (2026-11-07): `monitors::ArtPlanes` (`inc/monitors/art_planes.h`) with separate reflection (near y) and transmission (far y) planes and one `flux(side)`, which replaces `Reflectance_XZ` / `Transmittance_XZ` (the same integral). Fixes the shared-array bug (§7). The far-field point monitors read the transmission plane (forward) and the reflection plane (backward).
   - [x] Volume DFT (2026-11-07): `monitors::VolumeDft` (`inc/monitors/volume_dft.h`) owns the running E/H volume DFT, the normalised complex fields, the P_NL band of the exact multipoles and the DFG forward copy (kept across sweeps, as before). `allocate` keeps the per-array conditions; `release` frees all four monitor sets in one place (the old block missed the volume arrays of `ABS_SIGMA`-only runs). The 30 legacy names are references into it. Bit-identical: exact multipoles with SHG, `volume_dft` (3D E/H output), `sigma_absorption`, regression cases, hydro split, layout 222.
   - [x] Monitor list (2026-11-07): `monitors::Monitor` (`inc/monitors/monitor.h`) with `sample(SampleTime{e, h, step})`; `build_monitors` resolves the run's flags once per sweep (in `Driver::BuildSweepLists`, with the sources) into the volume DFT, P_NL DFT, flux boxes, display planes, ART planes, field trackers and custom DFT boxes/probes, in the old order. `DFT_Update` and the volume-DFT free functions are gone (`VolumeDft::sample_e/sample_h`); the field-tracker loop is one `add` per component. Bit-identical against the previous build: field trackers (linear, and the hydrodynamic split with the `_nl` trackers), custom DFT box + time probe, display planes, volume DFT / multipoles / absorption, periodic suite, regression cases, hydro split, layout 222, unit tests.
   - [x] Normalisation per monitor (2026-11-07): `Monitor::normalize()`; the driver calls it on the list where it called `DFT_Norm`. One `monitors::normalize_plane(normal, re, im)` (`inc/monitors/plane_dft.h`) replaces `DFTNormalizeXY/XZ/YZ` (three copies differing in the plane axes; missing components skipped); `VolumeDft::normalize(E/H)` and `normalize_pnl` replace `DFTNormalize` and `NLPolarizationDFTNormalize`; `DftPlanes`, `FluxBoxes` and `ArtPlanes` normalise their own planes under the conditions of `DFT_Norm`. `normalizeDFT.cc` (334 lines) is removed. Bit-identical with phase normalisation on and off: all display planes with H, volume DFT with H, exact multipoles with the P_NL band, cylinder SHG (forward/backward far-field planes, Q and efficiency files), periodic ART, hydro split; regression cases, layout 222, unit tests. Not covered by a run: the graphene plane.
   - [x] Analysis and output per monitor (2026-11-07): `Monitor::finish()`, called on the list after the normalisation with a barrier after each: the volume DFT runs the exact multipoles, the ART planes `ART_PBC`, a near-to-far monitor writes the scattering-box faces, the flux boxes the Q coefficients and ART-TFSF files, a field-point monitor the near/far-field point files (its flags moved out of the driver). The finish phase is collective, so the list is now built from run-wide flags only and the rank-local conditions (all-CPML ranks, tracking points) are tested inside `sample`/`normalize`. The always-true `dft_this_sweep` is removed. Bit-identical: the normalisation set (phase on/off), trackers and custom monitors, periodic suite, all sources, CPML-edge layouts with volume DFT + trackers + planes, regression cases, hydro split, layouts 222/411, unit tests.
   - [x] Output layout and plane maps (2026-11-07): `io::OutputLayout` (`inc/io/output_layout.h`: start, global, local, shift per axis) computed by `io::output_layout()` replaces the four driver-context arrays; the display-plane maps (`DFT_Plot`, or the volume E cut along the planes with field saves) moved from `FinishSweep` into the display-plane monitor's `finish`. `FinishSweep` is now the TTM pump file, the normalise/finish loop and the DFG bookkeeping. Bit-identical: planes with and without the volume save, the normalisation set, trackers, custom monitors, all sources, mesh dump, regression cases, hydro split, layouts 222/411, unit tests.
   - [x] End-of-run outputs (2026-11-07): `Monitor::write_final()`, called by `WriteOutputs` once per run (after the last sweep, after an interrupt, or on the mode-source dummy run): the 3D volume E/H files and the custom DFT boxes. `WriteOutputs` is now the loop plus the memory release. Bit-identical: volume saves with and without planes, custom monitors, mode-source dummy run with volume save, all sources, regression cases, hydro split, layout 222, unit tests.
   - [x] Field trackers (2026-11-07): `monitors::FieldTrackers` (`inc/monitors/field_trackers.h`) owns the points (local index and sphere-grid index) and `dft[set][E/H][component]`; `setup` replaces `GetMultipolePositions`, `sample` the per-component block, `write` the checkpoint files. 18 globals removed (no legacy references needed: nothing else read them). The DFG forward copy is `VolumeDft::keep_forward()`. Bit-identical: trackers linear and with the hydrodynamic harmonic set, trackers across CPML-edge layouts (x, y, z, 2x2x2) against one rank, DFG optimisation (two sweeps, volume save), custom monitors, regression cases, hydro split, unit tests.
   - [x] Time-domain probes and checkpoints (2026-11-07): `monitors::TimeProbes` (`inc/monitors/time_probes.h`) takes the `[probe.*]` points from the configuration, keeps this rank's cells and their E/H traces; it replaces `SetupTimeDomainProbes`, the `PROBES` table and its global (`apply_sim_config` no longer takes it), and the traces in the driver context. `Monitor::every_step()` marks monitors that record during the CW settling periods (the probes; the DFT monitors wait as before). `Monitor::checkpoint(step)` replaces the probe, tracker and Q/ART-TFSF blocks of `periodic_analysis`. Bit-identical: probes on 2x1x1 and 2x2x2, a CW run with settling (traces from step 1), probes with trackers, the normalisation set, custom monitors, all sources, regression cases, hydro split, unit tests.
   - [x] Custom monitor boxes (2026-11-07): `monitors::CustomBoxes` (`inc/monitors/custom_boxes.h`) owns the frequency and time boxes of the `[monitor.*]` sections (one `CustomBox` type; the `Monitor<T>` template, whose time instantiation carried unused DFT arrays, is gone with `Monitors.h`, `dcomplex`, `readMonitors.cc`, `setupMonitors.cc`, `writeMonitors.cc` and `UpdateCustomMonitorDFTs`); `build`, `write_masks`, `sample`, `write_dft`. The vectors and the flag left the driver context. The boxes are rebuilt per sweep from scratch (§7). Bit-identical: boxes on 2x2x2, 3x1x1, a rank edge on the CPML interface, two-sweep DFG, probes and trackers, regression cases, hydro split, unit tests.
   - [x] ART intermediate spectra (2026-11-07): `Monitor::after_step(step)`, called every step where `periodic_analysis` evaluated ART at steps N - kN/10; the ART monitor keeps that schedule (`IsArtCheckpoint` removed from the driver context). Bit-identical: periodic suite (all `ART_90pol_<n>` files), phase normalisation set, regression cases, unit tests.
   - [x] Near-to-far and exact multipoles in the module (2026-11-07): `monitors::write_near_to_far()` (`inc/monitors/near_to_far.h`) writes the scattering-box faces from `FluxBoxes` with one loop over the normal axis (replaces the three `Save_Surface_Near_to_Far_*` copies, `Write_Ini_File_NF` and `writeNearFar.cc`); `exactMultipoles.cc` moved to `src/monitors/exact_multipoles.cc` in `fdtd::monitors` (`sample_nonlinear_polarisation`, `write_exact_multipoles`). Pitfall found on the way: a function first named `exact_multipoles` in `fdtd::monitors` hid the `fdtd::state` flag of the same name inside the namespace, so `|| exact_multipoles` tested a function address (always true) and every run sampled an unallocated volume DFT (segfault in the display-plane runs). The compiler warned (-Waddress); the build output was filtered for errors only. Rule since: module functions are verbs (`write_*`, `sample_*`), and new modules are built with the warnings shown. Bit-identical: near-to-far faces and parameter file (2x1x1, 2x2x2, 1x3x1 with phase normalisation; datasets and attributes), exact multipoles with SHG (phase on/off), normalisation set, trackers, custom monitors, all sources, regression cases, hydro split, unit tests.
   - [x] ART and Q coefficients in the module (2026-11-07): `art.cc` -> `src/monitors/art.cc` (`write_art_spectra`, `integrate_absorption`), `qCoefficient.cc` -> `src/monitors/q_coefficients.cc` (`write_q_coefficients`, `write_art_tfsf`), in `fdtd::monitors` with headers; stale `Reflectance_CTW` / `Transmittance_CTW` declarations removed. Bit-identical: periodic suite, periodic Drude background with sigma absorption, normalisation set, `sigma_absorption`, regression cases, hydro split, unit tests.
   - [x] Field-point files (2026-11-07): `monitors::write_field_points()` (`inc/monitors/field_points.h`): one writer for the five copies of `updateFieldMonitors.cc` (removed), the five flag globals replaced by the conditions at the call. The near-field file takes the xz plane, else the yz plane; with neither DFT plane it is no longer written (it held uninitialised values, §7). Bit-identical: dipole with periodic ART (forward, backward, near), dipole with xz and with yz planes, cylinder SHG forward/backward, normalisation set, regression cases, unit tests.
   - [x] Plane DFT and plane maps (2026-11-07): `monitors::sample_plane(normal, pad, E/H, set, t, re, im)` (`inc/monitors/plane_dft.h`, with `normalize_plane`) replaces the twelve `*DiscreteFourierTransformUpdate{YZ,XZ,XY}[_NL]` entry points and `updateDFT.cc`; the display planes, flux boxes and ART planes call it with their own arrays (the function-pointer tables and their alias are gone). `DftPlanes::write_maps()` replaces `DFT_Plot` and `writeDFT.cc` (fifteen per-plane array lists). `src/monitors/` now holds every monitor: 14 files. Bit-identical: display planes with E/H and the harmonic planes (hydro split), the normalisation set, field points, periodic suite, near-to-far, all sources, regression cases, CPML-edge layouts with all planes vs one rank, layout 222, unit tests.
   - [x] Monitor legacy names removed (2026-11-07): the 312 `fdtd::state` references into the monitor objects are gone; their last users (the absorption integral, the exact multipoles, the plane maps) read the module objects. `integrate_absorption` is one sum over the components on the material (seven written-out cases before; same order of summation). `globals.h` no longer includes the monitor headers. The binary shrank from 2.6 to 1.7 MB (each `inline` reference had a dynamic initialiser in every translation unit); run time unchanged (100^3 Drude+2CP sphere). Bit-identical: normalisation set, periodic suite, periodic sigma absorption, volume DFT / multipoles / absorption, trackers and custom monitors, near-to-far, all sources, field points, custom boxes, regression cases, hydro split, unit tests.
   - [x] LWA multipoles (2026-11-13): a general analysis of the multipole monitor (`[multipoles] lwa`), see §7.
   - [x] CPML and auxiliary-grid legacy names removed (2026-11-07): the 70 CPML names (`psi_*`, `coeff_*_cpml`, `inv_stretch_factor_*`) and 44 of the 45 auxiliary-grid names (`ex_uni`, `hz_uni_pbc`, `psi_*_uni_*`, ...) are gone; their users (`update_h.cc`, `e_faces.cc`, `update_e.cc`, `hydrodynamics.cc`, `setupSimulation.cc`, `aux_grid.cc`, `tfsf.cc`) read the module objects through one reference cached per file (the accessors are out of line, the lesson of the `cpml()` slowdown). Kept: `mat_y_uni` (the 1D material map, read by the geometry, the hydrodynamic scans and the E scan). `fdtd::state` now holds 600 variables and 29 references (27 field components, `mat_y_uni`, `subdomain`). Binary 1.4 MB. Bit-identical: every material model as structure and background, polarisation 0/45/90 in all backgrounds, layered-medium waves, CPML-edge layouts in vacuum and Drude background vs one rank, layouts 222/411/141/114/311, periodic suite, normalisation set, all sources, regression cases, hydro split, unit tests; timing equal (vacuum and Drude background).
6. [~] **Free functions into their modules' namespaces;** the last `fdtd::state` aliases removed (done 2026-11-13).
   - [x] `grid` and `update` (2026-11-07): `SetFDTDlatticeData`, `SetAxisIndices`, `GetNeighboursRanks[NL]`, `GetProcessRank`, `GetProcessRankByCommunicator`, `PrintDecompositionLength` are `grid::set_lattice_data`, `axis_indices`, `neighbour_ranks[_nl]`, `process_rank`, `rank_in_communicator`, `print_decomposition` (declared in `inc/grid/subdomain.h`); `MultiSwitch` / `MultiSwitch2` are `grid::exchange_margins` / `exchange_margin` (`inc/grid/halo_exchange.h`); `ElectricFieldUpdateScan` / `MagneticFieldUpdateScan` are `update::update_e` / `update_h` (`inc/update/update.h`). Pitfall: a block-scope function declaration inside a function defined as `fdtd::update::update_e` declares `fdtd::update::ex_update_scan` (link error); such declarations go to file scope. Bit-identical: regression cases, hydro split, layouts 222/411, CPML-edge layouts in a Drude background vs one rank, periodic suite, every material model, unit tests.
   - [x] `materials` (2026-11-07): `MaterialsParameters`, `RefractiveIndex`, `ComplexPermittivity`, `ADE_Coefficients` are `materials::set_parameters`, `refractive_index`, `complex_permittivity`, `compute_ade_coefficients` (`inc/materials/materials.h`, by value instead of const references). Bit-identical: every material model, periodic suite, sigma absorption (vacuum and periodic Drude background), normalisation set, regression cases, hydro split, unit tests.
   - [x] `nonlinear` (2026-11-07): new module directory (`inc/materials/nonlinear/nonlinear.h`, `src/materials/nonlinear/{tensor,exact}.cc`, from `tensorNonlinear.cc`, `isotropicNonlinear.cc`, `inc/nonlinear.h`), namespace `fdtd::materials::nonlinear`: `Pxyz_NL_UpdateScan` / `Exyz_NL_UpdateScan` are `update_polarisation` / `update_field`, `NonlinearTerm` is `polarisation_term` (`polarisation` alone would hide a state variable), and the six `E{x,y,z}_Update_{SHG,THG}_DL` are `update_exact_shg/thg(a, ...)` (the component's field and its chi diagonal `[7a]`). Bit-identical: SHG/THG matrix (P_NL and field approach, staircased, exact multipoles), dielectric exact SHG and THG (chi verified active: a zero-chi run differs), regression cases, hydro split, normalisation set, unit tests.
   - [x] Signals and source setup into `sources` (2026-11-07): `signal.cc`, `pulseDFT.cc`, `modeSource.cc`, `readSources.cc` moved to `src/sources/{signal,pulse_dft,mode_source,plane_waves}.cc`; `Signal` / `dSignal` are `sources::signal_value` / `signal_derivative`, `PulseTransform` / `PulsePolarization` are `transform_pulse` / `select_polarisation`, the mode-source functions `setup_mode_source`, `write_mode_domain`, `read_mode_results`, `BuildPlaneWaves` is `build_plane_waves` (`inc/sources/signal.h`; `pulse`, `envelope`, `carrier` would hide state variables). The envelope/carrier helpers stay file-level until the signal becomes an object. Bit-identical: all sources, polarisation 0/45/90, layered waves, Gaussian beam / superposition, mode-source run, two-sweep DFG, regression cases, hydro split, unit tests.
   - [x] Geometry builders (2026-11-07): `basicGeometries.cc`, `legacyGeometries.cc`, `complexGeometries.cc`, `readInGeometry.cc` moved to `src/geometry/{basic,legacy,composite,read_in}.cc`; the 35 builders are `fdtd::geometry::create_*` / `read_geometries` / `create_structures` (`inc/geometry/builders.h`; `CreateRoatatedMonopole` is `create_rotated_monopole`). The set-arithmetic shapes keep their namespace `geom` (`inc/geometry.h`) until the builders are ported to them (§8). Bit-identical: the material maps of all 47 named shapes meshed on 2 ranks (`/tmp/meshall.sh`), regression cases, hydro split, unit tests.
   - [x] Set-up and I/O helpers (2026-11-07): `setupSimulation.cc` -> `src/app/setup.cc` and `memoryRelease.cc` -> `src/app/release.cc` (`app::setup_domain`, `initialise_run`, `allocate_fields`, `allocate_planes_and_material_maps`, `allocate_aux_source`, `allocate_scattering_and_absorption`, `release_run_memory`; `inc/app/setup.h`); `writeUtils.cc` / `textFile.cc` -> `src/io/` (`io::write_log`, `io::read_text_entry`; `inc/io/text.h`); `readInMagnitization.cc` -> `src/materials/magnetisation.cc` (`materials::read_magnetisation_field`). Bit-identical: regression cases, hydro split, every material model, all sources, custom boxes (two-sweep DFG), mode-source run, unit tests.
   - [x] Electron-gas and thermal functions (2026-11-07): `hydrodynamics.cc`, `specialFunctions.cc`, `writeMultipoles.cc` -> `src/materials/hydro/{hydrodynamics,electron_gas,lwa_multipoles}.cc` in `fdtd::materials::hydro` (`update_hydrodynamics`, `update_electron_density[_local]`, `update_p{x,y,z}[_2]`, `setup_eps_inf_{x,y,z}`, `update_nonlocal`, `write_lwa_multipoles`; `inc/materials/hydro/hydro.h`); `setupTemperature.cc`, `twoTemperature.cc`, `writeTemperature.cc` -> `src/materials/thermal/{setup,two_temperature,output}.cc` in `fdtd::materials::thermal` (`update_ttm`, `update_ttm_rk`, the temperature-dependent constants, `read_temperature`, `save_temperature`, ...; `inc/materials/thermal/thermal.h`; `load_temperature` would hide a state flag). Renames only: the per-component hydrodynamic copies and the TTM crash wait for their reviews. Bit-identical: hydro split (with and without the multipole sphere), every material model (nonlocal, GNOR, Drude-Lorentz), SHG/THG matrix, regression cases, unit tests; the TTM run fails the same way as before.
   - [x] `globals.h` now declares no functions: it is the shared state (`fdtd::state`) and the module includes. Step 6 continues with the state itself: the remaining variables move into the module that owns them (§2.4).
   - [x] State owned by one module moved out of `globals.h` into `inc/<module>/state.h` / `src/<module>/state.cc` (namespace `fdtd::<module>`): hydro (77), sources (40), monitors (34), geometry (18), materials (17), grid (12), thermal (6). `DipoleInput` moved to Sources.h. 396 declarations remain in `globals.h`; these are shared by several modules and need an owner chosen by hand.
   - [x] Layers for the shared state (2026-11-07, `docs/dev/StateLayers.md`): each remaining variable is owned by the layer that matches when it is written. Run configuration (fixed after configure) goes to the application layer; rank decomposition to `grid::Subdomain`; mesh products (material maps, model coefficients and state, chi tables, temperatures) to the module that builds them; run state (time, sweep, polarisation) to the driver context or the step arguments. Not one application struct for everything: per-rank, per-mesh and per-step values in it would recreate `fdtd::state`.
   - [x] `app::RunConfig` (2026-11-07, `inc/app/run_config.h`, `src/app/run_config.cc`): the 154 values fixed for the run, in nine groups (`time`, `grid`, `source`, `materials`, `nonlinear`, `thermal`, `geometry`, `analysis`, `output`). Read as `run_config()` (const, inline; `.cc` files cache `run_cfg` once); written only by the configure phases through `edit_run_config()` (`apply_sim_config`, `setup_domain`, `initialise_run`, `load_configuration`, `initialize_ranks`, `assign_output_paths`); the driver seals it after `initialize_ranks`, after which a write fails to compile and a late `edit_run_config()` throws. The classification was checked by the compiler: the uses were rewritten from the build errors, so every write outside the configure phases would have been an error (there were none). Bit-identical: all suites (regression cases, hydro split, monitors, custom boxes, probes, near-to-far, field points, normalisation, sources, periodic, mesh of all shapes, material models, polarisation, layered waves, SHG/THG matrix), unit tests; timing equal (100^3 Drude+2CP sphere and Drude background).
   - [x] `pids` retired (2026-11-07): `apply_sim_config` fills `RunConfig` from `SimConfig` directly, with no `INI_DATA_STRUCT` in between (873 `pids->` uses). Its 90 fields joined the `RunConfig` groups; where a field duplicated a `RunConfig` value (`WIDTH`, `LENGTH`, `pbc_x`, `save_output`, `ABS_SIGMA`, `pnl_approach`, ...) the two are one member, and the input codes that were not duplicates got names (`source.kind`, `analysis.dft_plane_mask`, `nonlinear.process`, `nonlinear.chi2_antenna`, ...). The per-rank cell counts of the decomposition (`nr_cells_x/y/z[]`) are `grid::cells_per_rank_x/y/z`. Bit-identical: all suites, unit tests.
   - [x] Shape flags (2026-11-07): the 46 booleans that each meant `antenna_type == shape::X` (34 in `RunConfig`, 12 in `fdtd::geometry`) are gone, with their 140-line set-up chain; the code tests `run_cfg.geometry.is(shape::sphere)`. Bit-identical: all suites, unit tests.
   - [x] Rank-decomposition copies (2026-11-07): 54 globals that duplicated `Subdomain` fields are gone; the readers use the field (TF/SF box, ART planes, flux boxes, animation planes) or a new `Subdomain` accessor (`nea_cpml_x()` ... from `which_type_of_process_*`, replacing a 60-line if-chain). The copies were conditional (0 when the feature was off); an instrumented run compared each copy with its field on the regression layouts: equal wherever the feature is on, and the readers run only then. Bit-identical: all suites, CPML-edge layouts (vacuum and Drude background) vs one rank, unit tests; timing equal (Drude background, sphere).
   - [x] Per-rank values into `Subdomain` (2026-11-07): `startindices`, `localsize`, `j_for_uni`, `has_plane_wave` are members (`output_start`, `output_size`, ...), set where they were set. Bit-identical: regression cases, hydro split, monitors, normalisation set, all sources, mesh of all shapes, unit tests.
   - [x] `rank_this`, `global_procs_size`, the halo datatypes and the output communicators into `grid` (2026-11-07): `grid::world_rank`, `grid::world_size`, ... Bit-identical: regression cases, hydro split, monitors, near-to-far, periodic suite, unit tests. The rank-decomposition layer (§2 of `StateLayers.md`) is done.
   - [x] Mesh products with a clear owner (2026-11-07): temperatures to `thermal`, chi tables / nonlinear flags and cell maps to `nonlinear` (new `inc/materials/nonlinear/state.h`), electron-gas arrays and constants to `hydro`; 61 variables, uses qualified from the build errors. Bit-identical: regression cases, hydro split, every material model (nonlocal, GNOR), SHG/THG matrix, unit tests.
   - [x] Last shared variables (2026-11-11): the material parameter tables, material maps, polarisations, poles and currents are `fdtd::materials` state; the per-sweep polarisation, the signal normalisation, the analysed frequencies and the superposition focus / box positions are `fdtd::sources`; `chx/chy/chz` are `fdtd::update`; the last hydro and monitor values went to their modules. The run state (time, step, sweep, step count, start clock, failure code, DFG flag, the two trace files) is `app::RunState` (`inc/app/run_state.h`), read as `fdtd::app::run_state().timestep`. `inc/globals.h` now declares no variables: `fdtd::state` holds only the 29 module references (field components, `mat_y_uni`, `subdomain`). Bit-identical: regression cases, hydro split, material models, all sources, monitors, periodic suite, SHG/THG matrix, custom boxes (two-sweep DFG), unit tests.
   - [x] `fdtd::state` removed (2026-11-13): the 29 references live in their owner modules (`fdtd::grid::subdomain`, `fdtd::fields::Ex` ..., `fdtd::sources::mat_y_uni`); the 74 `using namespace fdtd::state;` lines became per-file using-declarations of exactly the names each file uses (generated from the compiler's suggestions); the model headers use `fdtd::constants` only. `inc/app/globals.h` only gathers module headers now.
   - [x] Time as arguments (2026-11-13): no module reads the run's time or step any more. `StepTime` carries the step's time (`t`) for the mode source; `update_ttm(step)`, `save_temperature(step)`, `write_art_spectra(step)`, `write_lwa_multipoles(steps)` take it from the time loop; the ART monitor keeps the last step it saw for its final spectra (a run stopped by a signal named them one step short).
   - [x] `app/globals.h` removed (2026-11-13): every file that included it now includes the module headers it uses (found by compiling each file without it); `tools/globals_doc/` (the source data of its comments) removed.
   - [x] The dispersive state and the model parameters as objects (2026-11-13):
     `materials::DispersiveState` (`inc/materials/dispersive_state.h`; per E component P and its
     history, the critical-point poles, J, the hydrodynamic J history, the harmonic copies) and
     `materials::ModelParameters` (permittivities, conductivities, the E coefficients, Drude and pole
     constants). The kernel names (`Px`, `px_cp`, `jx`, `hydro::jx_prev`, `epre_inf`, ...) are
     references into them, so no kernel changed. The allocation, the memory report and the
     box-integrity check walk the object (`for_each_field`, `for_each_pole`) instead of their own
     lists. Bit-identical on every suite, plus a dispersive restart and the custom-box layouts; run
     time unchanged (within noise).
   - [x] (2026-11-13) The component view addresses the state through the component
     (`Component<A>::state()`, a `ComponentState&`); the nonlocal / hydrodynamic per-component arrays
     (pole sum, nonlocal P and its copy, grad div P and its history, the padded E copies, the
     electron density, the background permittivity) and the hydrodynamic parameters (nonlocal and
     GNOR coefficients, beta, densities, masses) are in the objects; `release_run_memory` and the
     restart file walk the object (`release()`, `for_each_saved`, same file order). Bit-identical,
     restarts of plain, Drude+CP and Lorentz runs included.
   - [x] The kernels use the members (2026-11-13): the 50 reference names of `materials/state.h`
     (`Px`, `px_cp`, `jx`, `epre_inf`, `ome_plasma`, `big_omega_cp`, ...) are gone; each function
     takes `mp = model_parameters()` (const unless it writes) or `ds = dispersive_state()` and
     reads `mp.epre_inf[mm]`, `ds.c[0].p` (about 650 uses in 31 files, rewritten by a token-level
     script; comments keep the formula names). Bit-identical: all suites. Hot loops timed against the
     commit before (production flags, 100^3, 1 rank, 300 steps, minimum of 3 to 6 alternating runs):
     dielectric, Drude, Drude+2CP, Lorentz, modified Lorentz backgrounds, SHG / THG field and
     polarisation updates, hydrodynamic split: all within the run-to-run noise (-5 % .. +2 %).
   - [ ] Next: the hydrodynamic model's own per-component names (`hydro::px_nl`, `grad_div_x`,
     `ex_2`, ...; `materials/hydro/state.h`) go when its x/y/z copies are folded (physics review).
7. [ ] Local items still open in §5-§8, taken up when their module is moved.
8. [x] **Layers (§2.5, done 2026-11-13: no core module reads the configuration or run state, no include goes against the order, `make check-layers` strict for all):** include hygiene (done 2026-11-13: `app/globals.h` removed, each file includes what it uses), the report `tools/check_layers.py` / `make check-layers` (done, counts in §2.5), `fdtd::legacy` removed, `grid::Lattice` / `GlobalPlanes`, `update::Spec`, `Cpml::background` (done; grid, fields, boundary and update read no configuration), then the core modules take typed parameters instead of `RunConfig`, bottom up.

---

## 4. Done: driver rework and legacy cleanup

- [x] INI configuration with command-line overrides; `mesh` / `check` / `simulate` commands; state-machine driver (Boost.SML, kept by decision).
- [x] Validated rank decomposition including 1 rank; `Subdomain` object per rank (119 `[MAXGLOBALPROCSSIZE]` tables and 47 copied globals removed); output communicators split on own membership.
- [x] carray ownership; `view3`/`view4` kernel parameters; HDF5 output.
- [x] No JSON input (geometry, plane waves, monitors, DFT boxes are INI sections; json-c dependency gone).
- [x] snake_case for 1371 globals and lattice fields; 67 duplicated values merged; integer codes replaced by names (`inc/kinds.h`).
- [x] Globals:
  - 42 unreferenced and 90 write-only globals removed, with the code that fed them (the hydrodynamic moment sums ran every step for nothing);
  - about 450 unused locals and the dead allocations removed;
  - about 1000 lines of commented-out code removed;
  - 130 single-file globals made file-local;
  - every remaining variable documented (`///<`, units, index order; source data in `tools/globals_doc/`, removed with `globals.h`);
  - all of them in `fdtd::state`. The namespace exposed two hidden collisions: the `DipoleSource` table vs the struct, and `count` vs `std::count`.
- [x] Comments: banner blocks removed, last Italian notes translated.
- [x] Naming decision: physics names keep their capitals (`Ex`, `Hy`, `DX`, `DT`, `OMEGA`); genuinely legacy names get renamed.
- [x] Geometry:
  - layers and set operators (`geometry.h`); the anapole;
  - sphere, cylinder, dipole, bowtie, triangle and the rounded family as set arithmetic. `CreateRoundMonopole` went from 673 lines to 12;
  - cross sections taken from the shape;
  - shape scope reduced to the keep-list: primitives, classic nanoantennas, metasurface/LIDAR. The photonic crystal, Brillouin and 28 other shapes are removed.

## 5. Done: kernels and performance

- [x] **Stretch scheme for the E scans** (`stretch_update.h`):
  - rows are split into stretches of one CPML region and one material;
  - dielectric and Drude+2CP stretches are branch-free;
  - the stretch lists are built once per mesh;
  - the CPML/TF-SF terms are a callback, so every local Drude+2CP cell takes the stretch path;
  - one copy of the Drude+2CP arithmetic (`Drude2CPPoint`).
- [x] **H scans:** one generic update per component, one `Psi` for all CPML faces of E and H (`updateH.cc` 795 -> 381 lines, `cmpl.cc` 1123 -> 402).
- [x] **Material models folded per component:** Drude+2CP, Drude, Lorentz, magnetic Lorentz, modified Lorentz, nonlinear Drude-Lorentz and the hydrodynamic per-cell update. The shared CPML/TF-SF terms are `CpmlTfsfTerms_X/Y/Z`.
- [x] **Nonlinear:**
  - one `NonlinearTerm` chi contraction;
  - one polarisation loop and one E-approach loop, over cell lists built once per mesh;
  - the dielectric exact SHG/THG is one quadratic and one cubic solver;
  - `tensorNonlinear.cc` 1395 -> about 200 lines.
- [x] **DFT planes:** one `plane_dft`, which gathers each plane once and then does a contiguous accumulate per frequency (`kernels::dft_accumulate`, optional DAMM backend).

| Case | Before | After |
|---|---|---|
| 120^3, 40 nm Ag sphere, 100 steps (E stretches) | 11.1 s | 7.2 s |
| 100^3 Ag sphere (stretch lists) | 7.30 s | 6.02 s |
| 100^3, Drude background (CPML terms in the stretch) | 17.9 s | 13.8 s |
| 100^3 SHG (nonlinear cell lists) | 7.23 s | 6.55 s |
| 64-frequency DFT run | 3.10 s | 2.38 s |

**DAMM findings** (kept for reference):
- composed DAMM operators are 2.5x slower than one fused row loop for the curl (memory-bound, five passes);
- the compiler already vectorises the row loop;
- DAMM bug: `fused_union` computes `rem_cols = M % tile_cols` where it should use `N`. Patch in `docs/dev/damm-fused_union-rem_cols.patch`.

**Still open** (most of these become part of §3 step 4):
- [x] The other models as policies, chosen once per stretch (§3 step 4, 2026-11-07).
- [x] Dispersive state limited to the material's bounding box (2026-11-07, §3 step 4); full per-stretch packing not needed for the memory.
- [ ] Hydrodynamic per-component scans (`SetupEPSINF*`, `P*UpdateScan(2)`, about 1700 lines). The X/Y/Z copies differ by 240-345 lines each, so they need a physics review before folding.
- [ ] Remaining strided surface loops: plane-wave / Gaussian-beam injection on XY faces, and the field-monitor `[m][...]` accumulators.
- [ ] Pitfall: `make clean` when switching `MARCH`. Mixed FMA and non-FMA objects give 1-ulp differences.

## 6. CPML review

- [x] **Dipole test** (`test/cpml/check.py`), 8 cells, probes 4 cells from the CPML: face -79.4 dB, corner -66.7 dB. Identical on 1 and 2 ranks and with a rank boundary inside the CPML.
- [x] **Plane-wave test** (`test/cpml/plane_wave.py`), normal incidence:

  | Background | CPML cells | Reflection |
  |---|---|---|
  | vacuum | 4 / 8 / 16 | -48 / -83 / -101 dB |
  | n = 1.5 | 8 | -79 dB |
  | n = 3.5 | 8 / 16 | -60 / -90 dB |

  In dense media the reflection follows the cells per wavelength, so the grading needs no medium correction. Dividing sigma by n (the old `N_AVG` lines) was tested and rejected: it would multiply sigma by about 6 in a metal background.
- [x] **Fixed:** the CFS alpha grading at E positions (integer `m / N`); a dipole source on a rank-boundary plane being dropped; a segfault on ranks that are all CPML along z (introduced with the stretch scans).
- [x] The three per-axis CPML range blocks are one `axis_cpml_ranges`, with the sentinels named.
- [x] **Sentinels:** `min`/`max` are true local indices except on ranks that are all CPML along the axis (empty range, see §3 step 3); `test/edge_layouts.sh` puts rank edges on every interface. Left: settle on half-open ranges at the call sites.
- [x] The `*_extra` values named `*_nea_cpml_end` / `*_far_cpml_start` (2026-11-13).
- [x] **Split `SetFDTDlatticeData`** (2026-11-13): eight named steps of `LatticeSetup` (§3 step 3).

## 7. Known bugs and open questions

| Area | Status |
|---|---|
| Two-temperature model | `Update_TTM` crashes (pre-existing). Reproduce with `[model] temperature = true`. |
| Split / code 200 | Some split-field runs abort with code 200. Make it a clear `fdtd check` error. |
| Hydrodynamics | Single-wavelength runs write empty `Qcoeff` tables (also upstream). |
| Complex build | Removed 2026-11-07 (it did not compile); complex and single-precision fields return as a template parameter (§2.2). |
| `grating` geometry | Dead code path: delete it or wire it to a shape name. |
| Layered-medium injector | Confirm 1 rank against 8 ranks bit for bit. |
| `WIDTH`/`LENGTH`/`LENGTH2` | One geometry rescales the global while the cross-section normalisation reads the original. The shapes moving to declared cross sections removes this. |
| Magnetic Lorentz background | `geometry.first/middle/last_medium = 400` segfaults (pre-existing; same in the build before the E-scan merge). As a structure it runs. |
| Changed 2026-11-19 | Materials by name, the entry above. Changes of results, none in the suites: (1) the ITO Drude materials take the plasma frequency from the carrier density with the solver's electron charge and mass, `E0`, `ME`; the table had its own literals, so the plasma frequency moves by 5e-10 relative. (2) The two-temperature plasma frequency starts from the heated material's own value; the table took, for films under 50 nm, slot 17, which is empty, so the offset was the negative of the temperature fit. (3) A material with a tensor is nonlinear in every shape and every part it fills; the table made the structure nonlinear for a list of shapes only, and the `bulk` mode made a cube of the layers nonlinear, which is gone: give the layer's material a tensor. (4) `[electrons] diffusion` was given to codes 100, 101, 300 and 304 only; it is a key of each material. (5) A dispersive material used only by a block, the shell or an anapole fill had no state allocated before; every slot has. (6) `[nonlinear] graphene_optimization` no longer overrides the indices of the dielectrics 1, 2 and 3. |
| Open 2026-11-19 | `refractive_index` and `complex_permittivity` of a Lorentz material differ: the first has 2 i omega gamma in the denominator, the second i omega gamma. Both kept as they were, `src/materials/permittivity.cc`; one of them is wrong. The hydrodynamic, Drude-Lorentz and magnetic Lorentz models have no permittivity formula there and read 0, so a medium of those models gives the sources and the Q normalisation an index of 0. |
| Fixed 2026-11-21 | Harmonic efficiency (`harmonic_efficiency*.txt`), a known error since the legacy code: (1) the input energy was the Parseval sum of the source DFT over the pump band only, so the pump spectrum outside the band was missed and the efficiency too high: the pump band held 52 % of a 20 fs pulse at 1.26 um over 1210 to 1310 nm, a factor 7.5 for the broadband pulse of `base_nl.ini` over 400 to 700 nm; (2) both band sums gave every point a full step, edges included, a rectangle rule biased by half a step at each edge, the energies 4 % apart between 21 and 81 points. Now the input fluence is n/eta0 times the sum of E_inc^2 dt over the source signal, and the harmonic energy is integrated with the trapezoid rule or, `spacing = gauss`, Gauss-Legendre. The efficiencies of earlier runs are not comparable. |
| Changed 2026-11-21 | HDF5 outputs stored x fastest: `[z][y][x]`, planes `[z][x]`, `[y][x]`, `[z][y]`; `numpy.transpose` gives the old order. Scripts that index by `[x][z]` change; the electron temperature input of the two-temperature model is read in the same order, so a file written before this change has to be transposed. |
| Changed 2026-11-21 | `cylinder_shell` is drawn along y, the plane-wave propagation axis; it was along z. Same keys, `core_length` and `shell_length` now along y. Results of cylinder_shell runs change; no suite case uses the shape. |
| Open 2026-11-20 | Nonlinear dielectric with off-diagonal tensor entries, unconfirmed from reading the code: a non-dispersive slot takes the exact per-cell solve only when its chi2 or chi3 is diagonal-only (`iso_chi_2`, `iso_chi_3`, set in `src/materials/slots.cc`); `nonlinear::update_field` and `update_polarisation` take only dispersive slots (`nonlinear_model`). A dielectric with any off-diagonal entry, an isotropic chi3 for example, then falls through to the linear branch of `src/update/update_e.cc` and its nonlinearity is dropped without a message. Confirm with a test; route it through the field update or reject it in `fdtd check`. The old name "ISO_CHI3, isotropic (has only diagonal nonzero elements)" means diagonal-only, not isotropic. |
| Open 2026-11-20 | `[nonlinear] update = polarization` (pre-existing, same in the baseline `Pxyz_NL_UpdateScan`): `P -= P_NL` on the stored ADE state every step acts as a force -P_NL / dt^2 in the second-order recursion, so the effective nonlinear polarisation is P_NL / ((Omega0^2 - Omega^2) dt^2). Drude has no restoring term: gain 1 / (Omega dt)^2, 4.0e4 at 1.26 um with dx = 2 nm (checked on the recursion), so an ITO chi3 of 2.3e-18 blows up near 2e6 V/m (cylinder_shell anapole run, NaN at step 10800; reproduced at 4 nm, also with the poles fixed and with a 100x weaker pulse; `update = field` stable). Modified Lorentz (`P_cp[0] -= P_NL / dt`) gets about 1 / (Omega0 dt)^2, 470 for the Si3N4 UV pole. Grows as 1/dt^2, so it does not converge. Fix: inject (Omega0 dt)^2 P_NL for Lorentz poles, reject Drude-type slots. Both updates take the same `polarisation_term` and a constant chi3: neither is dispersive. |
| Open 2026-11-20 | Run tooling from the same run: (1) no check for non-finite fields, the run went on 19000 steps after the NaN; stop at the first one in the trace. (2) A catalog file `<name>.ini` with sections of another material (`[material.silver_d2cp.pole.1]` in `ito_d2cp.ini`) loads without them and without a message; make it an error. (3) `fdtd memory` reports "dispersive state written outside its material box" on 4x4x4 (`examples/anapole.ini`): lazy mode leaves the shared zero row of a boxed array unfilled and glibc can return reused memory for it; zero that row always. False alarm, real runs zero every array. |
| Magnetic Lorentz outside the TF/SF box | Neither its update nor the E-face pass adds CPML psi or TF/SF terms for it (Drude and Lorentz get theirs from the E-face pass and `sources::tfsf_e`). Only usable inside the TF/SF box until it takes the boundary callback. |
| Bloch / constant transverse wavenumber | Removed 2026-11-07 (it needed complex fields and the config already rejected it): `[grid] bloch`, the `[bloch]` section, `PBC_CTW`, the two-sweep H/E injection, `Setup_PBC_CTW`, the Bloch phase factors and 20 globals; also the always-true `forward_pulse` and its backward-pulse branches. Plain periodic x/z boundaries are unchanged. To be rebuilt later with complex fields (§2.2). |
| `center_freq` | Holds a wavelength (m). Rename it when the sources move (§3 step 5). |
| Fixed 2026-11-12 | CW pulse signal (`cw_pulse`): the carrier period was read uninitialised in `signal_value`; it is 1 / fundamental frequency, as the envelope and the step count already assumed. |
| Changed 2026-11-12 | Output paths and file names (new tree, frame names), `Qcoeff` headers, one polarisation per run: a run with `polarization` other than 90 used to run twice (the second polarisation defaulted to 90) and now runs once. Data otherwise bit-identical (regression cases, hydro split, monitors, custom boxes, SHG/THG matrix, sources, periodic, near-to-far, field points, probes, material models). |
| Fixed 2026-11-12 | Circular polarisation: `[source] circular` was reset during domain setup, so every run was linear. Now `[plane_wave] polarization_axis = rhcp / lhcp`, with the quarter-period shift taken from the carrier of the signal. The [source] keys moved to per-type sections (config change, data bit-identical). |
| Open 2026-11-13 | Hydrodynamic density update, SHG source terms and the critical-point switch (old Advanced Model codes 13, 14, 15): the keys `update_density`, `shg_source`, `critical_point` were read by nothing, so these codes ran like code 3. The keys are removed; restore the physics with the Phase 2 model review if it is wanted. |
| Changed 2026-11-13 | LWA multipoles are a general analysis: `[multipoles] lwa = true` writes `Multipoles_lwa_double_<pol>.txt` (and `_full_current_`, `Multipoles_lwa_crossx_*`) from the same volume current as the exact set, as its kr -> 0 limit (Alaee et al. Table 1), for every material model. The hydrodynamic LWA writer it replaces (`hydro/lwa_multipoles.cc`, 60 time-series arrays) had no code filling the arrays and its only call needed the always-false `e_dens_ft`; it is removed. The exact outputs are unchanged. Check: a 5 nm sphere, LWA within 2e-4 (700 nm) to 6e-4 (400 nm) of the exact moments, the (kr)^2 scaling. |
| Fixed 2026-11-13 | The derivative of the `cw_pulse` signal (used by `[wave_superposition] derivative = true`) read an uninitialised period; it now uses the fundamental, as the signal itself does. |
| Fixed 2026-11-13 | Release of the dispersive state: `release_run_memory` (called once, at the end of the run or after a mesh-only run) reset only part of the state, and six `if (*_nl != NULL)` lines left from deleted statements chained the temperature and energy-density resets onto the harmonic fields being allocated. It now releases the whole state object and the temperature arrays unconditionally. No output changes (nothing runs after it). |
| Open 2026-11-13 | `gaussian_alt` adds a second Gaussian whose carrier frequency is never set (0): the second pulse is the bare envelope, not a second colour. Decide what the signal should be. |
| Fixed 2026-11-07 | Drude+2CP coefficients: the sums of b2, b3, b4 over the critical-point poles that enter c0..c3 were initialised once before the loop over all material codes, so every Drude+2CP material after the first slot (gold, 100) also summed the poles of the materials before it. Silver (101) and every later Drude+2CP material had wrong E-update coefficients; gold was right. Now per material. Gold runs unchanged (regression cases bit-identical); a silver sphere (36^3, 400 steps) changes by up to 7.6 % in the scattering and 1.2 % in the absorption cross section. |
| Fixed 2026-11-07 | Near-field point file (`NEAR_FIELD_<pol>.txt`): with the animation planes on but neither the xz nor the yz DFT plane, the writer printed uninitialised amplitudes and phases; the file is now written only when one of the planes exists. |
| Fixed 2026-11-07 | Custom monitor boxes, multi-sweep runs (DFG optimisation, double polarisation): each sweep appended the boxes again, so sweep 1 held every box twice (each written twice to the same file, one communicator per copy never freed), and DFT boxes on ranks outside a box were sized from unset counts. The boxes are now rebuilt from scratch per sweep (communicators freed) and DFT arrays exist only where the rank holds a part. Outputs unchanged. |
| Fixed 2026-11-07 | Periodic ART: the near (reflection) and far (transmission) planes accumulated into one set of arrays, so on a rank holding both (every run with one rank along y) R and T were integrated over their sum and came out as T = -R (A = 1 exactly). The planes are separate now: a 2x1x1 vacuum run gives the same R and T as 1x2x1 (where the planes were already on different ranks, and which is unchanged). The far-field point monitors (forward/backward) read the matching plane. |
| Fixed 2026-11-07 | Dipole source, last step: the precomputed table had `N` rows but the loop runs steps 1..N, so step N read one row past it (0 on this machine, so the last step injected nothing); now N + 1 entries, which changes the final step of dipole runs. 1D simulation output: the last column of `1D_SIMULATION/*.txt` interpolated one past the 1D arrays (garbage, different per run); the profile ends one point earlier. |
| Fixed 2026-11-07 | Gaussian beam NaN: `[gaussian_beam] focal_factor` and `waist_factor` defaulted to 0, giving a zero focal length and waist and 0/0 amplitudes; they (and `theta_max`) are now required > 0 and documented. With values set the beam runs, bit-identical to the build before the superposition fold (layouts 211, 222). |
| Fixed 2026-11-06 | χ table collision (dispersive vs dielectric slot: an Au/Ag antenna in vacuum had χ = 0, upstream 92ea34d, report sent); hydrodynamic Ex CPML coefficient (`c_cp_1` -> `c_cp_0`, was NaN in the CPML); magnetic Lorentz Ey magnetisation (`magnetization_y` -> `_x`). |

## 8. Smaller open items

- [x] **XDMF for ParaView** (2026-11-07): every HDF5 output gets a `.xmf` companion (`[output] xdmf`, default on), written by `h5io::File` when the file closes from the datasets it wrote: 2D and 3D datasets on uniform grids with the run's cell size, datasets of one extent on one grid, the axis mapping in a comment. Complex data is now two real datasets `<name>_re` / `<name>_im` instead of the HDF5 compound (ParaView's reader cannot read compounds; decided with the user). Checked: VTK's `vtkXdmfReader` (ParaView's reader) loads the plane, volume, near-to-far, mask and custom-box files with the right dimensions and spacing, and its arrays equal the HDF5 data; every value of the split datasets equals the compound values of the previous build (normalisation set, all sources, near-to-far, custom boxes, probes, regression cases); `test_h5io` checks the split and the companion. `test/regression/run_tests.py` reads both layouts. User page: `docs/tutorial/Outputs.md`.
- [x] **`definitions.h` removed** (2026-11-07): `using namespace std;` gone from every file (the code qualifies `std::`; unqualified `abs` made `std::abs`); `REAL`/`REAL2` are `double`, `CMPLX` is `std::complex<double>`, the MPI type macros are `MPI_DOUBLE`/`MPI_DOUBLE_COMPLEX`, `RWCAST`/`ROUNDNEARESTINTG` written out. Constants: `inc/constants.h` (`fdtd::constants`); material codes and their layout: `inc/materials/codes.h`; the legacy parameter block, probe block and failure codes: `inc/legacy_parameters.h` (`fdtd::legacy`, transitional). The complex-field build (`PRECISION=complex`, did not compile) and its 35 `#ifdef DOUBLECOMPLEX` blocks are removed; precision and complex fields come back as a template parameter (§2.2). Bit-identical, mesh included.

- [x] **Compiler warnings** (2026-11-13): a clean `-Wall -Wextra` build (solver, tests, sigcfg) has none. Fixed: two assignments used as conditions (decomposition, intended; parenthesised as the others), the TF/SF-face rank indices read uninitialised when no rank holds a face (now "no limit"), `catch (std::bad_alloc)` by value, the aggregate init of `LatticeSetup` (member initialisers), signed/unsigned compares; removed unused and set-but-unused locals (about 150 assignments in the hydrodynamic per-component copies, left from removed code paths); unused builder parameters are commented in the signature (`/*l1_l*/`), which keeps the record of which keys a composite ignores. Also `make clean` now removes the objects of the module subfolders (it left `src/materials/*/*.o`, so a clean build reused stale objects).
- [x] **Nonlinear tensors in the INI** (2026-11-13): the chi2 / chi3 tensors of `[nonlinear] chi_structure / chi_first / chi_middle / chi_last` are `[chi2.<name>]` / `[chi3.<name>]` sections (non-zero components by name: `d11` ... `d36` in m/V, `xxxx` ... `zyzy` in m^2/V^2) or built-in presets with the values of the former `INIDEF/chi2|chi3/*.txt` files (`config/chi_tensors.h`; every preset checked value by value against its file before the files were removed). The tensor-file reader and the files are gone; `[paths] input_dir` holds data files only (mode profiles, electron temperature). Intentional changes: an unset `chi_*` name is the zero tensor (all four were required; at least one must be set); a name with neither a section nor a preset is a validation error (`fdtd check` reports it; it was a file-open error at run time). Bit-identical: nonlinear suites (864 / 172 files) with the test tensors moved from files into the INI.
- [x] **`Doxyfile` untracked** (2026-11-13): `make doxygen` generates it (`doxygen -g`) and sets the project lines; the ones the tracked file carried are Makefile variables (`DOXYGEN_INPUT`, `DOXYGEN_MAINPAGE`, `DOXYGEN_IMAGE_PATH` for the shape figures) plus `RECURSIVE = YES` for the module folders; output in `docs/output`. `make cleanall` removes the generated file. Not run here (no doxygen on this machine); the overrides were checked on a copy of the old file.
- [x] **Function names to snake_case** (2026-11-13): every function and method (178 names, about 960 uses in 75 files; code, comments and docs). The list came from a libclang inventory of all declarations, so only names that are functions and nothing else were renamed by token; 147 took the plain conversion (`LoadSimConfig` -> `load_sim_config`, `ExUpdateScan` -> `ex_update_scan`, `SaveAsBinary2D` -> `save_as_binary_2d`). The 31 whose plain form is a keyword, a type, or an existing variable or function got a name by what they do: the row set operations `Union` / `Intersection` / `Difference` -> `unite` / `intersect` / `subtract`, `QNormalisationArea()` -> `compute_q_normalisation_area()`, `sm::Run` -> `run_machine`, `DftPlaneMask` -> `dft_plane_mask_bits`, `geom::Anapole` -> `anapole_parts`, the three `Psi` kernels -> `psi_update`, file-local helpers (`Plane` -> `correct_plane` / `exchange_plane`, `At` -> `face_value`, `Make` -> `make_shape`, ...); the `Driver` phases are `decompose`, `prepare_run`, `finish_run`, `write_outputs`, ... Types stay PascalCase. Bit-identical: all suites; unit tests; no warnings.
- [x] **The remaining legacy variable names** (2026-11-13): the upper-case and mixed-case variables are snake_case (711 uses in 46 files, token-level script). `RunConfig` members: `DT`/`DX`/`DY`/`DZ` -> `dt`/`dx`/`dy`/`dz`, `WIDTH`/`LENGTH`/`LENGTH2`/`THICKNESS`/`MIDDLE_LAYER`/`ANTENNA_SHIFT`/`SHELLMATERIAL` -> `width`/`length`/`length_2`/`thickness`/`middle_layer`/`antenna_shift`/`shell_material`, `TFSF`, `DIPOLE_SOURCE`, `PULSEWIDTH` -> `pulse_width`, `NONLOCAL`, `GNOR`, `Hydrodynamics`, `CARRIER_DENSITY(_BACK)`, `Pump`, `Probe`, `ABS_SIGMA`, `ART`, `MULTIPOLES`, `PLOT_XY/XZ/YZ`, `TIME_PLOT_3D` lower-cased; the booleans `nonlinear.SHG`/`THG` are `second_harmonic`/`third_harmonic` (`shg`/`thg` are the mode codes). `ModelParameters::NR_CP` -> `nr_cp`. Module variables got names that cannot be shadowed by a local: `OMEGA` -> `analysis_omega`, `MAXFREQUENCY` -> `max_source_frequency`, `DFT_XY/XZ/YZ` -> `dft_plane_xy/xz/yz`, `Abs_sigma` -> `abs_sigma_spectrum`, `Area` -> `q_normalisation_area`, `GAP` -> `antenna_gap`, `SHG_bulk`/`THG_bulk`, `CHI2` -> `chi_2`, `PumpProbeDelay` -> `pump_probe_delay`, `JPWSTARTUNI` -> `jpw_start_uni`, `JPWSTARTUNI_` -> `jpw_source_uni`, `IPW1_`/`JPW1_`/`KPW1_` -> `ipw_near_h`/`jpw_near_h`/`kpw_near_h`. Constants keep upper case (`EP0`, `ROOTPROCESS`, `NRMAXMATERIALS`). Left: upper-case locals inside the hydrodynamic and structure code (`TOTND1`, `POINTRADIUS`, ...), with their owners. Bit-identical: all suites.
- [ ] **Doxygen headers for the functions.** They come with the module headers, not as a separate pass.
- [x] **Commented-out code removed** (2026-11-13): about 1800 lines in 53 files (`//` lines ending in `;`, braces, `if`/`for` heads, preprocessor lines; `/* */` blocks that were mostly code; runs of `//` lines with no prose left). Doc comments (`///`, `/**`) and prose comments are kept. The machine code differs only in string addresses and `__LINE__` constants; every suite is bit-identical.
- [x] **Multipoles as a post-processing choice** (2026-11-13): exact and long-wavelength both enabled from `[multipoles]` for any run, not tied to the hydrodynamic branch.
- [x] **Memory before the run: `fdtd memory`** (2026-11-14). One process, no `mpirun`, counts the arrays every rank of the configured layout allocates, without allocating them. It is not a table of costs per cell: `Driver::estimate_memory` plays each rank in turn through the solver's own decomposition, per-rank initialisation, allocations and mesh, and reads a tally kept by the one allocator (`palign`, `fdtd::core::memory`, bytes per group); in lazy mode `palign` maps address space, so only what the set-up writes costs memory. A 1000³, 64-rank configuration (124 GiB of arrays) takes 11 s and 280 MiB. The same tally is printed in the geometry report of `fdtd mesh` (per group, largest rank and total) and as `arrays: peak` at the end of a run. Checked: `fdtd memory` equals the peak of the real run to the printed digit in 53 configurations of the working suites (every material model, sources, monitors, periodic, SHG / THG, hydrodynamic split, two-temperature, magnetic Lorentz, blocks, volume DFT, uneven rank layouts), and `make regression` compares it with `fdtd mesh` on three cases. Not counted: the process overhead (a constant, 60 MiB per rank; 40 to 65 measured) and output buffers at the end of a run. Found on the way: `hydro::matind_n` was allocated (4 bytes per interior cell) and never used, removed; `release_run_memory()` did not release the magnetisation, `nd_z_prev`, `electron_temp_tf` and the multipole trackers (leaks between sweeps), fixed, and an array that outlives a rank is now an error of `fdtd memory`; the signal files were written through unchecked `fopen` results. The reports are in GiB (1024³ bytes). Bit-identical: all suites.
- [x] **Documentation pass** (2026-11-14): every page of the user guide describes the INI configuration and the current output tree. Rewritten: `README.md`, the home page, Running a Simulation, The Configuration File, Source Types, Plane Wave Superposition, Testing, the documentation guide; new: `models/Materials.md` (material codes), `models/NonlinearTensors.md`, `geometry/Blocks.md`, `dev/Architecture.md`, `dev/Migration.md` (the old-format key and code tables, formerly `docs/config_keys.md`); usage sections of the model, Gaussian-beam, dispersion and layered-wave pages converted; every INI example passes `fdtd check`. Removed: the onboarding documents, the cluster pages, `.gitlab-ci.yml`, `makepphxl`, the JSON geometry page, `tools/snake_rename.py` and its map. `test/compare_runs.py` is the bit-for-bit comparison of two output trees used for every refactoring; the `check.py` scripts take the output directory. `examples/anapole.ini` and `dimer_shg.ini` had too few time steps to pass `fdtd check` (60 and 40000; now 6000 and 60000). The three GitHub workflows (cmake build, `doxygen` without a Doxyfile, the GitLab mirror) were removed as well. Noted, not changed: several material slots have historical labels (`models/Materials.md`).
- [x] **Style pass and build** (2026-11-14): no banner comments, no asides in parentheses in the comments of `src`, `inc`, `test`, `tools` and the INI files (about 1500 lines, by script, formulas and the shape sketches kept); private members are `_member` (63 names); the `=====` rule under the geometry report title is gone. `sigcfg` is now the command `fdtd signal`, so the Makefile has one target. The Makefile only sets variables and adds `check-tests`, `regression`, `check-layers`; `rules.mk` gained `MARCH ?=`, dependency files (`-MMD -MP`, included and cleaned) and nothing project specific. HDF5: `pkg-config $(HDF5_PKG)`, default `hdf5`; `io/h5io.h` stops the compilation when the HDF5 found has no MPI support (`H5_HAVE_PARALLEL`). Bit-identical: all suites against the previous commit, `make regression` 10/10.
- [x] **Anapole keys and declarations** (2026-11-14): `[anapole]` defaults are the published 200 nm design (Hassan et al. 2025), which is the size the lobe fit refers to; `corner_radius` 11 nm, read off the paper's figures, the fit gave none. Keys renamed from the mirrored COMSOL drawing to what they are on the antenna: `gap_length`, `gap_width` (full sizes, were `neck_half_length`, `neck_half_height`), `lobe_length`, `lobe_height`, `lobe_exponent`, `lobe_offset` (were `lobe_a`, `lobe_b`, `lobe_p`, `lobe_x0`). Fills: `aperture_material` (lobes and gap) and `gap_material`; the dimer is removed. The regression case sets the 400 nm variant explicitly and matches its stored reference, so the rename is bit-identical. `examples/anapole.ini` is the published case with the gap monitor and the variants, `examples/anapole_ito.ini` the aperture filled with ITO (code 153). Check runs of the default shape: dip at 800 nm with 2 nm cells and 880 nm with 4 nm, against 735 nm in the paper; see `docs/geometry/Anapole.md`. Declarations: one per line (504 split, found with libclang), one space between type and name. Bit-identical: all suites, `make regression` 10/10.
- [x] **Materials as definitions in slots** (2026-11-19): the built-in table of material codes is gone, `materials/parameters.cc`, `materials/codes.h`, `materials/inputs.h`, `config/chi_tensors.*`. A material is a definition, `[material.<name>]` with a `model` and its keys, poles in `[material.<name>.pole.<i>]`, `extends` for variants; tensors are `[chi2|chi3.<name>]`. The solver reads them from the run file or the catalog, text files in `materials/` and the directories of `[materials] catalog` and `$FDTD_MATERIALS`, root only before the broadcast; `resolved_config.ini` carries every definition used. A run loads its materials into slots, vacuum first; the material maps, the parameter, coefficient and state tables are indexed by slot, and the model is `kind(slot)` instead of a code range. Nonlinearity is a property of the material: a slot whose tensor is not zero is nonlinear wherever it is drawn, `[nonlinear] shg / thg` are booleans, the per-medium tensors and the nonlinear code offset are gone. Config keys removed: `[materials] n<i>`, `carrier_density*`, `stretch_x`, `[custom_material]`, `[nonlinear] chi_*`, `[electrons] diffusion`; `tools/migrate_materials.py` converts a file. The catalog holds every entry of the former table, checked value by value against a dump of the old table; `docs/models/Materials.md`. Bit-identical on all suites with the configurations converted by name; the material maps and masks differ by the renumbering of codes to slots only, checked map by map, and the regression references of the two maps are regenerated. Intentional changes, §7.
- [x] **Report names and metadata** (2026-11-21): the text reports have plain names without the polarisation or the precision, `cross_sections.txt`, `art_tfsf_<n>.txt`, `art_<n>.txt`, `harmonic_efficiency*.txt`, `multipoles_<method>*.txt`, `source_signal.txt`, `source_spectrum.txt`, `aux_grid_<field>.txt`, `far_field_*.txt`, `near_field.txt`, `run.log`, and in `frames/` `probe_<x>_<y>_<z>.txt` and `tracker_<n>.txt`; the full table is in `docs/tutorial/Outputs.md`. Each opens with `io::open_report`: comment lines with the version, precision, config file and overrides, source and polarisation, spectrum, grid and time step, then the column header, also as a comment, so `numpy.loadtxt(path)` reads the numbers. The header of `source_spectrum.txt` names its six columns, `source_signal.txt` has one; the first lines of the probe and tracker files, the material slots and the cell, are comments. `sources::pol` is gone; `write_log` writes one `run.log` and closes it. Numbers identical on all suites, file by file through the rename map; `make regression` with the references renamed, the unit tests, `check-layers`, no warnings.
- [x] **HDF5 on the solver's axes** (2026-11-21): every HDF5 file stores x fastest, volumes `[z][y][x]`, planes `[b][a]` (`XZ` as `[z][x]`), so a plane is an image with its first axis across and ParaView's X, Y, Z are the solver's. `h5io` takes blocks in the solver order and reorders through a dense buffer in 32-wide tiles (`for_each_reordered`), on write and read; the memory selection of the old path is gone. The `.xmf` grids are all 3D, a plane one point thick along its normal at its new `position` attribute (frames, material and DFT planes; the near-to-far faces at 0). Checked: every dataset of all suites equals the old one transposed (7,020, planes and volumes), text reports unchanged; `vtkXdmfReader` puts every point of a volume and of each plane on its solver coordinates on an anisotropic 40 x 30 x 50 grid; the frame, material and DFT planes of one run coincide; `test_h5io` checks the stored order on 1 to 8 ranks; `make regression` with the references transposed; the shape figures unchanged.
- [x] **Either spacing with a harmonic band** (2026-11-21): the input and harmonic energies of the efficiency reports sum each point times the frequency width it stands for, `sources::fundamental_point_width` and `monitors::harmonic_point_width` in place of one step: the step for a band uniform in frequency, as before; half the gaps to the neighbours for a band uniform in wavelength, an edge point the gap to its one neighbour. `[spectrum] spacing = wavelength` with a harmonic band is no longer rejected. Bit-identical on all suites. THG on `base_nl.ini`, 400 to 700 nm: the scattering efficiency with wavelength spacing is 4.2136e-9 and 4.2115e-9 at 21 and 81 points against 4.2111e-9 and 4.2109e-9 with frequency spacing.
- [x] **Harmonic efficiency: time-domain input, quadrature, Gauss spacing** (2026-11-21): the input fluence is n/eta0 times the sum of E_inc(t)^2 dt over the source signal, n averaged over the sampled pump spectrum; the harmonic energy uses `config::spectrum_weights`, trapezoid for uniform spacings, Gauss-Legendre for the new `spacing = gauss`, whose nodes `spectrum_omegas` places. Checked: on `base_nl.ini` with 120 Gauss nodes over 60 nm to 100 um the Parseval sum of the source DFT is 0.9945 of the time-domain integral, the rest outside that band; THG harmonic energy 6.667561e-44 J with 3, 5, 9 and 17 Gauss nodes, 6.667374e-44 with 5 and 6.667554e-44 with 21 uniform points; all suites identical except the 24 `harmonic_efficiency*.txt`; unit tests for the nodes, symmetric, exact to degree 2n - 1, and the weights summing to the band.
- [x] **Band estimates in `fdtd signal`** (2026-11-21): with a `[spectrum] harmonic_band`, `band:` lines report before the run the spectral FWHM of the pump and of its k-th power, the share of each inside its band, the error of the solver's band rule against a converged integral, a band centred on the carrier that holds `-c` (default 0.999) of the harmonic, the noise efficiency (the pump's own energy in the harmonic band, DFT over the configured steps at the solver's nodes and weights, over the input fluence) and the rounding floor; `-e` gives the margin of an expected efficiency. New `core/quadrature` (Gauss-Legendre, moved from `config`, and panels) and `sources::SampledSignal`; `test_spectrum` checks them against the closed forms of a Gaussian. The report now places the points with the solver's C0. Cylinder anapole, 20 fs at 1.26 um over 1210 to 1310 nm: pump 0.524 and harmonic 0.783 in band, the analytic values 0.524 and 0.783. Solver unchanged: `make regression` 10/10. See `docs/tutorial/SignalSettings.md`.
- [x] **`fdtd signal` spectrum table** (2026-11-21): the amplitudes are the DFT over the configured steps, as the monitors take it, instead of over the signal up to the 1e-3 threshold, whose cut raised the harmonic-band amplitudes to 2e-5 of the peak (2e-13 over 36000 steps for the cylinder anapole); `WEAK` marks pump points only, and harmonic-band points, which the solver does not normalise, are marked `harmonic`.
- [x] **Console output and reports, grep-able and concise** (2026-11-21): every console line starts with a key; data are tables, a key whose first line is the comma-separated header and whose next lines are the rows, so `grep '^key:'` is a CSV file: `fdtd signal` (`signal`, `spectrum`, `harmonic`, `floor`, `band`, `noise`, `config`), `fdtd memory` and the geometry report (`layout`, `memory`, `bytes`, `job`, `grid`, `time`, `resolution`, `peak`, `census`), `fdtd check` (`; check:`), and the solver's log (`start`, `domain`, `subdomain`, `planes`, `flux_y`, `n_avg`, `nonlinear`, `qcoeff` with a `step` column and one header for the intermediate and final passes, `q_area`, `q_norm`, `faces`, `art_area`, `anapole`, `mode`, `centroid`, `beam`, `trace`, `arrays`). Messages are short and without parentheses, warnings `warning:`, errors on stderr `error:`. Help texts are usage and options only. Reports carry one metadata line, program and configuration file, then snake_case column names with the unit last; memory groups are snake_case. Key descriptions of `fdtd template` and validation messages lose their parentheses. `test/regression/run_tests.py` reads the `memory:` and `bytes:` tables; reference headers renamed. Checked: `make regression` 10/10; the named suites against the build before, 1949 files, data identical.
- [x] **Upstream modules restored** (2026-11-19): `extern/` holds carray, upstream commit 678c610, and the INI parser, `config.h`/`config.cc` as received, unmodified; `exceptions.h` was not received and is reconstructed. The solver's additions moved into derived modules: `core::array` over `carray` (empty state, shape, boxed storage, memory charge, `field3`/`field4`/`vector1`; aliases `array1`..`array4`) and `io::Ini` over `cfg::Config` (parsing from text, case-insensitive names, `#` comments, boolean words, signed integers, sections in file order). `validation_error` moved to `config/errors.h`. The upstream patches that would remove each workaround are listed in `extern/README.md`. Memory reports identical to the byte; bit-identical on all suites, `make regression` 10/10.
- [x] **Named parameters per antenna** (2026-11-14): the last five shapes have a section, so every shape does and `[geometry]` holds no dimensions (`width`, `length`, `thickness`, `length_2`, `gap`, `arm_length`, `arm_width`, `rotation`, `corner_radius` are no longer keys; they remain the internal members the builders read). `[l_antenna] width, thickness, length_x, length_z`; `[lateral_dipole] width, thickness, gap, left_length, right_length, overlap`; `[round_butterfly] width, thickness, gap, length_z, length_x, corner_radius, rotation`; `[round_tripole] width, thickness, center_length, top_length, bottom_length, top_gap, bottom_gap, corner_radius, rotation`; `[x_shape] thickness, fixed_length, fixed_width, rotated_length, rotated_width, corner_radius, angle`. The names were chosen from the builders and checked on the mesh (the figures in `docs/geometry/Shapes.md` place each dimension line from the key value). Intentional change: `rotation = 0` now means bars along z for the tripole (builder angle 90) and bars along x and z for the butterfly (builder angle 45), counter-clockwise positive, as for the round dipoles; `tools/migrate_config_keys.py` converts old files, including the angle. Meshes identical to the previous build with the old keys (the five shapes, two rotations; the other 27 shapes unchanged).
- [ ] **Port the remaining kept composites** in `complexGeometries.cc` / `legacyGeometries.cc` to the geometry layers.
- [x] **Anapole regression case** (2026-11-13): `test/regression/Anapole`, the anapole nanoantenna as a 30 nm silver patch (code 101) in air on a glass substrate, plane wave from the air side, 4 nm cells, 8 ranks, 1200 steps (a regression check, not a converged spectrum). References: Q coefficients, ART-TFSF, the XZ material map and the E DFT plane at 1200 nm.
- [x] **Shipped regression references regenerated** (2026-11-13): `make regression` had not been run since the references were written on 2026-11-03 and 5 of the 6 cases failed. Causes: (1) the output tree and frame names changed on 2026-11-12, so the reference files were at the old paths (moved: `ART` -> `reports`, `DOMAIN` -> `materials`, `MOVIE/Efield_tscaled=<n>_XY.h5` -> `frames/Efield_XY_<n>.h5`); (2) the values moved by up to 1e-6 relative in T and R (ThinFilm, Mie) with the CFS alpha grading fix of 2026-11-06 (commit 72a85bb; found by bisecting the ThinFilm ART file, the build flags were ruled out); (3) the `Qcoeff` files have a header line since 2026-11-12; (4) `EmptyDomain_Gaussian` runs a pulse longer than the run and now needs `allow_truncated_signal = true`. The references were regenerated from the current build and pass on a second run at the 1e-10 tolerance. `make regression` is part of the checks from here on.
- [x] **Pre-/post-processing scripts** (2026-11-13, decided by the user): `preprocessing/`, `postprocessing/` and `SimulationSetup/` are deleted, including the near-to-far transform (`make post`) and `fdtd_io.py`; rebuild what turns out to be needed against the HDF5 outputs.

## 9. Review notes from the user (2026-11-12)

Requirements from a review of the branch. Ordered inside each phase from cheap to large; `[ ]` = to do.

### Phase 1

**Outputs**
- [x] `Qcoeff*` files have column headers (2026-11-12).
- [x] Movie frames are `frames/Efield_<plane>_<00000>.h5`: the sequence number last, zero-padded (2026-11-12).
- [x] Output tree (2026-11-12): `reports/`, `fields/`, `frames/`, `materials/` (the material index maps), `state/` (application state dumps); the eleven old directories and `OPT/` are gone, the 1D profiles are `reports/1D_*.txt`, and a directory the run leaves empty is removed at the end. The scripts in `postprocessing/` still use the old names (§8). Was: Target (decided 2026-11-12): `reports/` (the text tables, named by report), `fields/` (frequency-domain, DFT), `frames/` (time-domain movie frames), `state/`. `OPT/` goes (it does nothing now); `1D_SIMULATION/` only when the 1D simulation writes output.

**Logs**
- [x] Grep-able log (2026-11-12): `io::log(key, format, ...)` (`inc/io/text.h`) writes one record per line as `key: text` (`key[rank]: text` from a rank other than 0); blank lines and rules are dropped. All 154 `printf` sites go through it with their module as the key (`grid`, `setup`, `run`, `sources`, `monitors`, `materials`, `hydro`, `thermal`, `geometry`); tables have their own keys: `qcoeff:` / `qcoeff_nl:` (the Q spectra at each checkpoint), `n_avg:`, `trace:` (the field trace). `grep '^qcoeff:' log` extracts the table. Left: the wording of the older messages, and a step number on the `qcoeff` checkpoints.

**Signal and time steps**
- [x] Guard (2026-11-12): `sources::signal_fits` samples the signal over the run; the run stops with a message when the signal is still above 1e-3 of its peak in the last 2 % of the steps (or is zero throughout). Steady signals are exempt; `[source] allow_truncated_signal = true` overrides (short test runs).
- [x] `sigcfg` (2026-11-12, `tools/sigcfg.cc`, `make sigcfg`; `docs/tutorial/SignalSettings.md`): builds the signal as the solver does and reports its length, the minimum and recommended steps, the amplitude at the analysed wavelengths and the config line, as key-prefixed lines. For this the signal selection and timing moved out of `initialise_run` into `sources::configure_signal()` and the time step into `app::set_time_step()` (both without MPI). `fdtd check` runs the same test. The guard in `simulate` stays as a backstop: it runs first in the run preparation, before the field allocation, and costs two signal evaluations per step.

**Config**
- [x] Source and signal sections (2026-11-12): `[source] type` and `signal` select; the keys of a type are in the section of the same name (`[plane_wave]`, `[wave_superposition]`, `[gaussian_beam]`, `[dipole.<name>]`) and the keys of a signal in `[signal.<name>]` (`max_frequency` for the raised cosines and the rectangle; `center_wavelength`, `width`, `chirp*` for the Gaussians; `settle_periods`, `periods` for the CW signals). Setting a key in the section of a type or signal that the run does not use is an error (the section must equal its defaults, so the template and `resolved_config.ini` still reload). The waves of a superposition are `[wave.<name>]` (were `[plane_wave.<name>]`). `[source]` keeps `peak_field`, `allow_truncated_signal`, `aux_steps`, `aux_length_factor` (was `wave_factor`: it is the length of the 1D grid, not an amplitude). `force_aux` is `[wave_superposition] layered_simulation`. Removed: `circular`, `circular_delay`. `tools/migrate_config_keys.py` rewrites old files.
- [x] Electron and temperature models (2026-11-13): `[model]` is split into `[electrons] model` (+ `diffusion`, `fermi_velocity`, for the non-local models only), `[hydrodynamic]` (`pressure`, `convection`, `magnetic`, `split_nonlinear`), `[temperature] model`, `[two_temperature]` (`save_frames`, `rk_stages`) and `[pump_probe]` (`delay`), with the same rule: a key set in the section of a model the run does not use is an error. `[anapole]` is checked the same way against `[geometry] shape`. Removed, because nothing read them: `update_density`, `shg_source`, `critical_point` (§7). `tools/migrate_config_keys.py` rewrites old `[source]` and `[model]` sections.
- [x] Geometry keys cleaned (2026-11-13): `[nonlinear] dfg_frequency` / `dfg_step` are `[geometry] rotation` [deg] and `corner_radius` [m], which is what the shape builders read them as (`run_cfg.geometry.rotation_deg`, `geometry::corner_radius`; `meta_mosaic` reads the radius as its lattice step). Removed, read by nothing: `[geometry] height`, `ald_thickness`, `afm_size_x/z`, `[materials] stretch_z`, `sheet_spacing`. `[nonlinear] chi_*` set with `shg` and `thg` off is an error.
- [x] Shape sections, first group (2026-11-13; decided: keep the primitives and the nano-antennas). 21 shapes have a section named after them with keys named for what they are (`[dipole] width`, `thickness`, `gap`, `top_end`, `bottom_end`): the spheres (`sphere`, `sphere_shell`, `two_spheres`, `sphere_dimer`, `sphere_trimer`, `sphere_quadrumer`), `cylinder`, `cylinder_3d`, `dipole`, `bowtie`, `triangle`, `rounded_triangle`, `rounded_monopole`, `rounded_dipole`, `rounded_dipole_offset_planes`, `round_dipole`, `u_shape`, `h_shape`, `v_shape`, `l_shape`, `two_films_spacer` (and `anapole`, which had one). One table (`shape_sections()`, `src/config/sim_config.cc`) maps each key to the pool member its builder reads (`resolved_geometry`), so the builders are untouched and the meshes are identical. Setting a `[geometry]` dimension with such a shape, or a key in another shape's section, is an error. `[round_dipole] rotation` counts from the z axis (the builder counts from x).
  Figures: `docs/geometry/make_shape_figures.py` meshes every shape with `fdtd mesh` and draws the keys on the mesh as dimension lines placed from the key values (`docs/geometry/shapes/<shape>.png`), and writes `docs/geometry/Shapes.md` (figure and key table per shape, from `fdtd template`). The builder comments in `src/geometry/basic.cc` reference the figures (`@image html`) where they had ASCII drawings.
- [x] Shape sections, second group (2026-11-13): `cylinder_shell` is new (a `cylinder` core of `radius` inside a shell of `shell_thickness` and `[geometry] shell_material`); `hole` (`radius`, `thickness` of the film), `one_side_round_dipole` (`top_end`, `bottom_end` to the centres of the end discs) and `round_hybrid` (`arm_length`; the builder's formula `(length_2 + arm_width + width - gap) / 2` is inverted by the section's `finish` hook) have sections and figures; the `anapole` has its figure (with a close-up of the neck). The three rotated dipoles count `rotation` from the z axis.
- [x] 2026-11-13, decided by the user: `cylinder` and `cylinder_shell` are finite (`length` along z; 0 runs through the whole domain, which is how an infinite cylinder is configured), `cylinder_shell` is two coaxial cylinders with independent sizes (`radius_core`, `core_length`, `radius_shell`, `shell_length`), the shell drawn first and the core over it, so the shell can close the ends or leave them bare; `ellipse` has a builder (`create_ellipse`: an elliptic disc in the x-z plane, `radius_x`, `radius_z`, `thickness`, `rotation`; it was in the shape list without one); the shapes that did not draw what they are named for are removed: `z_shape` (filled most of the domain), `patch` (a film with slits; ignored its width and length), `bullseye` (one hole, no rings), with their builders and their branches in the background builder, the decomposition and the field points. `monitors::single_photon` was only set by the bullseye and is now always false.
- [x] Mode source rebuilt (2026-11-13): `[source] type = mode` with a `[mode]` section reads MPB's field files of one band (`efield`, `hfield`, `unit_length`, `effective_index`), interpolates them onto the Yee positions of a plane normal to x, y or z (`propagation_axis`, `position`, `mpb_x_axis`, `mpb_y_axis`, `offset_*`) and injects them as a one-way TF/SF plane (`src/sources/mode_source.cc`, `docs/sources/ModeSource.md`). Measured -46 dB upstream (MPB 20 nm, FDTD 40 nm), the same split along and across the axis and with the plane on a rank boundary. Replaces the two-pass flow (pre-run, `Mode.py`, text profiles) and the `mode` signal. Fixed on the way: the full material maps (PML included) that the mode-source update and CPML read were never filled (every cell used material 0), and the H correction was evaluated half a step late (-28 dB upstream). `[block.<name>]` sections (boxes of one material, file order, optional `count` / `pitch_*`) draw layered and periodic structures before the shape. Examples: `examples/mode_source/` (strip waveguide), `examples/grating_antenna/` (dipole on a Si3N4 grating coupler), each with its MPB input.
- [ ] Mode source: restart support; dispersive materials on the injection plane; a multi-frequency profile for broadband pulses.
- [ ] Other propagation axes and `propagation_angle` under one common frame: the TF/SF injection, the 1D grid and the ART planes are along y.
- [x] Shape sections for the remaining composite antennas: done 2026-11-14 (§8, named parameters per antenna).
  - `l_antenna`: `width` both arms; `length` the arm along x; `length_2` the arm along z; `thickness`. No rotation.
  - `round_butterfly`: `width`, `length_2`, `gap`, `thickness`, `rotation`, `corner_radius` (which also grows the extent along y); `length`, `arm_*` unused.
  - `round_tripole`: three collinear arms; `width`, `length`, `gap`, `arm_length`, `arm_width`, `thickness`, `rotation`, `corner_radius`; runs out of a 96 nm domain at 12 to 30 nm values.
  - `lateral_dipole`: two arms with diagonal cuts; `width`, `length`, `length_2`, `gap`, `arm_length`, `thickness`.
  - `x_shape`: `width`, `length`, `length_2`, `gap`, `thickness`, `rotation`, `corner_radius` (also along y).
- [x] `sigcfg` reports the pump's power at 2f and 3f for the configured steps (the floor a harmonic must rise above) and how it moves with the step count (2026-11-13).
- [x] Shapes that are not kept removed (2026-11-13): `meta_*` (metasurface lattices), `meta_cavity` (and its Ex/Ez trace files), `lidar_*`, with `[lidar]`, `[materials] damping_*` / `extra_layer_*` (read only by the LIDAR shapes), their branches in the plane placement, the background builder (the oxide / Ti / extra layers of the 1D grid), the material parameters (the Ti adhesion layer, the damping factors) and the dispersive-pole selection, and the builders only they used (`create_bolt_diag1_round`, `create_bolt_diag2_round`; `create_dbr` was unused). `tools/migrate_config_keys.py` drops the keys. The shape-code table in `docs/dev/Migration.md` is generated from `inc/geometry/shapes.h`.
- [x] `[custom_material]` is checked (2026-11-13): its keys are an error unless the structure, a layer, the shell, a `[block.*]` or an anapole part uses its code (302, the user hydrodynamic metal).

**File layout**
- [x] No file outside a module (2026-11-13): `src/` holds only `main.cc`, `inc/` only directories. Moves, with no code change:
  - `inc/core/`: `carray.h`, `views.h`, `aligned_memory.h`, `kernels.h`, `constants.h` (module-independent)
  - `config`: `ini.{h,cc}` (was `config`), `sim_config.{h,cc}`, `apply.{h,cc}`, `exceptions.h`
  - `app`: `driver.{h,cc}`, `driver_context.h`, `time_loop.cc`, `machine.h` (was `app_machine.h`), `globals.h` (removed 2026-11-13; the transitional `fdtd::state` references), `legacy_parameters.h`
  - `io`: `h5io.{h,cc}`, `binary.h` (was `io.h`)
  - `geometry`: `csg.{h,cc}` (was `geometry`), `anapole.{h,cc}`, `render.{h,cc}`
  - `sources`: `types.h` (was `Sources.h`), `gaussian_beam.cc` (was `setupSources.cc`), `one_dim_simulator.{h,cc}`, `plane_wave_one_dim_simulator.{h,cc}`
  - `materials`: `parameters.cc` (was `materials.cc`), `coefficients.cc` (was `adeCoefficients.cc`), `material_field.h`
  Removed: `inc/legacy.h` (an include of three headers), `src/globals.cc` (empty), the unused `DipoleSource` struct. The material-index dump was written out twice in the driver (run and mesh-only path) and is one function now.
- [x] The config layer is in `fdtd::config` (2026-11-13): the INI reader and its exceptions (was `cfg`), `SimConfig`, the section types, `load_sim_config` / `validate_sim_config` / `to_ini`, the chi tensors and `apply_sim_config` (were in `fdtd`).
- [ ] Left as they are: `OneDimensionalSimulator`, `save_as_binary_2d/3D` and `build_gaussian_plane_waves` keep their names and the global namespace; `ensure_directory` is only used inside `write_utils.cc`. The two 1D simulators (`sources/aux_grid` for the TF/SF plane wave, `one_dim_simulator` for superposed waves through layers) solve the same 1D problem with separate code: candidates to merge when the propagation axes are generalised.


### Phase 2

- [ ] Nonlinear processes and models reviewed with the material models (the current code makes material-agnostic assumptions); real time-domain wave mixing (two-colour source, DFG/SFG inside the solver) as a later item.
- [ ] Varin's nonlinear Lorentz drive on the Drude-Lorentz register (codes 350-399), as the SHG / THG path for dispersive materials. The register is Drude (the hydrodynamic or two-temperature electron gas) + one linear Lorentz pole kept in the critical-point arrays; it already integrates the left side of P'' + G P' + W^2 P = W^2 eps0 (chi1 E + chi2 E^2 + chi3 E^3) (Varin et al., CPC 222, 70 (2018); OME 9, 771 (2019)). Add the chi2 E^2 and chi3 E^3 drive with E at step n (explicit), from the chi2/ / chi3/ tensors or scalars, and compare with the nl_code tensor path on the dimer SHG example. The legacy chi_nl / NONLINLorentz1/2 coefficients were computed but never used (models/drude_lorentz.h).
- [ ] Materials: integrate the user's C++ shared library that fits oscillators to material data (Matadata), to replace the material-to-update-model mapping; Matadata may need anisotropic models. Review the shifters (the code-range offsets per model) and the material-table bloat that comes with them.
