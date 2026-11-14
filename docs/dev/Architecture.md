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

| Module | Owns |
|---|---|
| `core` | `carray` and views, aligned memory, constants, limits |
| `io` | `h5io` (parallel HDF5 with XDMF companions), text tables, logging |
| `grid` | `Lattice` (cells, PML, cell sizes, periodic faces), `Subdomain` (a rank's part of the grid), the decomposition, halo exchange, the CPML coordinate `Stretching`, the global plot planes |
| `fields` | E and H arrays and their halos |
| `materials` | material codes and the material map; the dispersion models (`models/`: Drude, Drude + critical points, Lorentz, modified Lorentz, Drude-Lorentz, magnetic Lorentz) with their parameters (`ModelParameters`) and state (`DispersiveState`); and three parts that act on a material: `nonlinear/` (χ² and χ³), `hydro/` (nonlocal and hydrodynamic electrons), `thermal/` (two-temperature model) |
| `boundary` | the CPML |
| `update` | the E and H update scans |
| `sources` | the time `Signal`; the sources as `Source` objects: TF/SF plane wave with its 1D auxiliary grid, plane-wave superposition and Gaussian beam with their 1D layered simulations, dipoles, the mode source |
| `monitors` | the monitors as `Monitor` objects: flux boxes and Q coefficients, ART planes, DFT planes, volume DFT, custom boxes, time probes, near-to-far faces, multipoles |
| `geometry` | the `Scene` to draw into; shapes as set operations on primitives (`csg.h`, `render.h`, the anapole) and the point-wise builders of the other shapes; blocks |
| `config` | the INI reader, `SimConfig` (one struct per section), validation, the template, the χ tensors, and `apply_sim_config`, which translates `SimConfig` into the application's `RunConfig` |
| `app` | `Driver` (the phases of a run) and its state machine, `RunConfig`, `RunState`, the set-up that hands each module its input, the time loop, state files for restart, signal handling |

## How a run is assembled

The application translates the configuration into one typed input per module. A module never looks
at the configuration itself.

| Input | For | Set in |
|---|---|---|
| `grid::Lattice`, `grid::GlobalPlanes`, `grid::Stretching` | `grid`, and everything above it | `Driver::decompose` |
| `update::Spec` | `update` | the set-up of the run |
| `materials::Inputs`, `materials::ModelSpec` | `materials` and its parts | the set-up of the run, after the time step is known |
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
| read, override and validate the configuration; broadcast it | `load_configuration` |
| decompose the domain over the ranks | `decompose`, `initialize_ranks` |
| build the signal, allocate, build the geometry, the sources and the monitors | `prepare_run` |
| `fdtd mesh` stops here | `dump_geometry` |
| continue an interrupted run | `load_state` |
| the time loop: H update, E update, sources, monitors | `run_time_loop` |
| normalise the DFTs and run the analysis | `finish_run` |
| write the bulk outputs | `write_outputs` |
| an interrupted loop writes its state instead | `save_state` |

## Memory accounting

Every array is allocated by `palign()` (`inc/core/aligned_memory.h`), which keeps a tally of the
bytes held per group (`fdtd::core::memory`); a `Scope` names the group of the allocations made while
it is alive. Three things read the tally: the geometry report of `fdtd mesh`, the `arrays: peak`
line at the end of a run, and `fdtd memory`. The last one plays every rank of the layout in one
process (`Driver::estimate_memory`): it sets the rank, runs the decomposition, the per-rank
initialisation, the allocations and the mesh, reads the tally, and releases everything. In its lazy
mode `palign()` maps address space instead of allocating, so only what the set-up writes (the
material maps of one rank at a time) costs memory. An array that survives the release between two
ranks is an error, which also keeps `release_run_memory()` complete.

## Conventions

- Functions, variables and configuration keys are snake_case; types are PascalCase
  ([Style Guide](StyleGuide.md)).
- A module's shared variables are in `inc/<module>/state.h` ([Shared State](StateLayers.md)).
- A refactoring leaves the results bit-identical; a change that moves results is recorded in
  `Roadmap.md` §7 with its size ([Testing](Tests.md)).
