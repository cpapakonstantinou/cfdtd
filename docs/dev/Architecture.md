# Architecture

The solver is a set of modules, one directory each under `inc/<module>/` (headers) and
`src/<module>/` (sources), in namespace `fdtd::<module>`. This page is the map; `Roadmap.md` §2 holds
the design reasoning and the record of how the code got here.

## Layers

| Layer | Modules | Role |
|---|---|---|
| 0 primitives | `core`, `io` | arrays and views, constants, kernels' building blocks; HDF5 (`h5io`), text and binary writers |
| 1 core library | `grid`, `fields`, `materials`, `boundary`, `update`, `sources`, `monitors`, `geometry` | the FDTD method |
| 2 application | `config`, `app`, `src/main.cc`, `tools/` | reads the INI file, assembles a run from the core library and drives it |

Two rules, checked by `make check-layers` (`tools/check_layers.py`):

1. A core module does not include application headers and does not read the application's run
   configuration or run state.
2. Core modules include each other only upward in this order:
   `grid < fields < materials < boundary < update < {sources, monitors, geometry}`.

## Modules

Upstream modules are kept unmodified in `extern/`: `carray` and the INI parser `cfg::Config`.
The solver extends them by derivation, `core::array` and `io::Ini`; `extern/README.md` lists
what the derived classes work around and the upstream patches that would remove it.

| Module | Owns |
|---|---|
| `core` | `array`, the solver's arrays derived from `carray`, the memory tally, views, constants, limits |
| `io` | `Ini`, the INI dialect on `cfg::Config`; `h5io` (parallel HDF5 with XDMF companions), text tables, logging |
| `grid` | `Lattice` (cells, PML, cell sizes, periodic faces), `Subdomain` (a rank's part of the grid), the decomposition, halo exchange, the CPML coordinate `Stretching`, the global plot planes |
| `fields` | E and H arrays and their halos |
| `materials` | material definitions, `Definition` read from INI sections; the slot table of the run's materials, `load_slots`, and the material maps of slot numbers; the dispersion models (`models/`: Drude, Drude + critical points, Lorentz, modified Lorentz, Drude-Lorentz, magnetic Lorentz) with their parameters (`ModelParameters`) and state (`DispersiveState`); and three parts that act on a material: `nonlinear/` (χ² and χ³), `hydro/` (nonlocal and hydrodynamic electrons), `thermal/` (two-temperature model) |
| `boundary` | the CPML |
| `update` | the E and H update scans |
| `sources` | the time `Signal`; the sources as `Source` objects: TF/SF plane wave with its 1D auxiliary grid, plane-wave superposition and Gaussian beam with their 1D layered simulations, dipoles, the mode source |
| `monitors` | the monitors as `Monitor` objects: flux boxes and Q coefficients, ART planes, DFT planes, volume DFT, custom boxes, time probes, near-to-far faces, multipoles |
| `geometry` | the `Scene` to draw into; shapes as set operations on primitives (`csg.h`, `render.h`, the anapole) and the point-wise builders of the other shapes; blocks |
| `config` | `SimConfig` (one struct per section), validation, the template, the material catalog, `with_catalog`, and `apply_sim_config`, which translates `SimConfig` into the application's `RunConfig` |
| `app` | `Driver` (the phases of a run) and its state machine, `RunConfig`, `RunState`, the set-up that hands each module its input, the time loop, state files for restart, signal handling |

## How a run is assembled

The application translates the configuration into one typed input per module. A module never looks
at the configuration itself.

| Input | For | Set in |
|---|---|---|
| `grid::Lattice`, `grid::GlobalPlanes`, `grid::Stretching` | `grid`, and everything above it | `Driver::decompose` |
| `update::Spec` | `update` | the set-up of the run |
| the `materials::Definition`s, `materials::LoadSpec`, `materials::ModelSpec` | `materials` and its parts | the set-up of the run: the slots first, the model spec after the time step |
| `geometry::Scene` | `geometry` | before the structure is built |
| `sources::Spec`, the `Signal` | `sources` | end of the run's initialisation |
| `monitors::Spec` | `monitors` | `Driver::configure_run` |

Where one core module needs something from a module above it in the order, the application passes
an interface instead of an include: the CPML receives the incident-field terms of the TF/SF source
as `boundary::IncidentETerms`.

## The phases of a run

`src/main.cc` runs the state machine of `inc/app/machine.h` over the `Driver`:

| Phase | `Driver` method |
|---|---|
| read the configuration and the catalog files it names, root only; broadcast; override and validate | `load_configuration` |
| decompose the domain over the ranks | `decompose`, `initialize_ranks` |
| build the signal, allocate, build the geometry, the sources and the monitors | `prepare_run` |
| `fdtd mesh` stops here | `dump_geometry` |
| continue an interrupted run | `load_state` |
| the time loop: H update, E update, sources, monitors | `run_time_loop` |
| normalise the DFTs and run the analysis | `finish_run` |
| write the bulk outputs | `write_outputs` |
| an interrupted loop writes its state instead | `save_state` |

## Materials

A material is data, `materials::Definition`: a model, `Kind`, and the model's parameters, read from
`[material.<name>]` and its `[material.<name>.pole.<i>]` sections by `materials::definition`, with
`extends` resolved. Each model lists its keys once, the key tables of `src/materials/definition.cc`;
reading, checking and writing a definition go through those tables. The configuration finds every
material the run names in the run file or in the catalog, `config::with_catalog`, root only, and
the broadcast text carries the catalog files.

The run loads its definitions into slots, `materials::load_slots`: slot 0 is vacuum, then the
structure, the three media, the shell, the anapole fills and the blocks, each material once. The
material maps hold slot numbers; the parameter, coefficient and state tables are indexed by slot,
and `kind(s)` is the model of slot s. A nonlinear material is a slot whose tensor is not zero; its
cells are nonlinear wherever it is drawn.

## Memory accounting

Every array charges the blocks it allocates when it is created (`inc/core/array.h`) and returns
the charge when its last copy releases the memory. The tally keeps the bytes held per group
(`fdtd::core::memory`, `inc/core/memory.h`); a `Scope` names the group of the charges made while it
is alive. Three things read the tally: the geometry report of `fdtd mesh`, the `arrays:` table
at the end of a run, and `fdtd memory`. The last one plays every rank of the layout in one process
(`Driver::estimate_memory`): it sets the rank, runs the decomposition, the per-rank initialisation,
the allocations and the mesh, reads the tally, and releases everything. In its lazy mode large blocks
come from fresh pages and are not zero-filled, so only what the set-up writes, the material maps of
one rank at a time, costs memory. An array that survives the release between two ranks is an error,
which also keeps `release_run_memory()` complete.

## Conventions

- Functions, variables and configuration keys are snake_case; types are PascalCase
  ([Style Guide](StyleGuide.md)).
- A module's shared variables are in `inc/<module>/state.h` ([Shared State](StateLayers.md)).
- A refactoring leaves the results bit-identical; a change that moves results is recorded in
  `Roadmap.md` §7 with its size ([Testing](Tests.md)).
