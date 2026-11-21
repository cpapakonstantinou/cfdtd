# Rework plan: profile, inventory, materials

Written 2026-11-19 for review, revised the same day after your answers. Step 1 of the order
agreed, carray and config restored as upstream modules in `extern/`, is done: a3b96a1, 623fdc3.
This document is step 2. Order as decided: the materials first, then module ownership.

Done: materials as definitions in slots, Roadmap §8. The catalog is a directory of text files,
`materials/`, outside the binary.

## Decisions

- Accuracy before speed: the DFT samples every step. No decimation.
- Report names as proposed in the table below, with the run's metadata as comment lines above the
  column header.
- matadata integration is deferred until you have worked on matadata; its dependencies stay out
  of the solver for now.
- The goal of the profile is the shape of the source, not the size of the binary: call chains that
  show where the script repeats itself or passes work through layers that add nothing.
- The materials table is redesigned: each material its own definition, loaded into slots.

## Measurements

- Case: the anapole at 10 nm, 50³ interior cells, 16-cell CPML on every side, 81 spectrum points,
  silver as Drude + 2CP, one rank. `/tmp/profrun/run.ini` is `examples/anapole.ini` on one rank.
- Cost: callgrind on `OPT="-O3 -g" MARCH=-march=x86-64-v3`, 150 steps, and the wall clock.
- Call chains: callgrind on `OPT="-O1 -g -fno-inline"`, 4 steps, so every source function appears
  as a call; `/home/claude/work/cg_graph.py` prints the tree and the forwarding functions.

Cost, for the record:

| part | share |
|---|---|
| flux-box DFT, memory bound | 65 % of the wall clock at 81 points |
| E update with CPML psi | 47 % of instructions |
| H update | 22 % of instructions |
| Q coefficients, normalisation | 6 %, 10 checkpoints in 150 steps |

The loops that matter are vectorised. The bit-identical speed items that remain are small and go
with the edits that touch them: 1/kappa along z as a per-k table, `src/update/update_e.cc:275`;
the two `MPI_Barrier` per step, `src/app/time_loop.cc:285,287`; one packed reduce in place of 12
per field set in `write_q_coefficients`.

## What the call chains show

**Work reaches its kernel through layers that only forward.** From the time loop to the work:

    Driver::update_electric > update::update_e > ey_update_scan > e_scan<1>
      > StretchList::for_row > row lambda > cell lambda > Model::point

    Driver::accumulate_monitors > FluxBoxMonitor::sample > FluxBoxes::sample > lambda
      > sample_plane > plane_dft > kernels::dft_accumulate

    Driver::periodic_analysis > FluxBoxMonitor::checkpoint > write_q_coefficients > face_flux

    FluxBoxMonitor::normalize > FluxBoxes::normalize > normalize_plane

`ex/ey/ez_update_scan`, `hx/hy/hz_update_scan`, the ten `Monitor` subclasses and the module
functions behind them each add a named function and no decision. The `Monitor` class is the
interface; the singletons behind it are the implementation; the wrapper is the cost of having
both.

**The application layer knows every module's lifecycle.** A monitor is named in five to nine
files: its own header and source, `monitors/monitor.cc`, `app/setup.cc` to allocate,
`app/state_file.cc` to save and restore, `app/release.cc` to free, `app/driver.cc` to order, and
the config. Flux boxes: 9 files; field trackers: 7; time probes, custom boxes: 5. Adding a feature
means editing all of them, which is why the code grows faster than the function.

**State is global, so the lifecycle has to be spelled out.** 151 `extern` variables in the module
state headers, 38 in `monitors/state.h`, 32 files marked `transitional`. Because the state is not
owned by anything, allocation, release, restart and the memory report each walk it by name.

**The materials table is the same pattern at its worst.** `ModelParameters` is a struct of arrays
indexed by code, `[600]` and `[500]`, the model is chosen by the code's range, `codes.h`, every
model has `coeff_table<M>[500]`, the state has `pc[500]` pole arrays per component, and
`set_parameters` is 1 086 lines: every material ever used as assignments into those arrays, and
the chi2 and chi3 checks written out four times, once per medium position.

## Design: modules own their state

The Linux driver model is a good guide here. A driver is a set of operations, `probe`, `remove`,
`suspend`, `resume`; the core holds a list of bound devices and calls the operations; resources
taken with `devm_` are released when the device goes, without the driver writing `remove`. In C++
terms:

| kernel | here |
|---|---|
| `struct device_driver` ops | a concept: `configure`, `allocate`, `sample`, `checkpoint`, `archive` |
| probe | the constructor, from the module's config section |
| `devm_` managed resources, remove | RAII: `core::array` members, `memory::Charge`; the destructor is the release |
| suspend, resume | `archive(Archive&)`, one function for save and load, as now |
| the bus's device list | the driver's tuple of modules, or a vector of monitors |

- Each module's state becomes one struct owned by the module, as `fields::Fields` and
  `grid::Stretching` already are. The externs go, and with them the `transitional` markers.
- A monitor is a class that satisfies the concept; the driver holds
  `std::vector<std::unique_ptr<Monitor>>` built from the config, and calls each operation in one
  loop. No wrapper classes, no per-monitor lines in setup, release or state_file.
- The fields, CPML, sources and material slots are always present, so the driver holds them by
  value in a `std::tuple` and `std::apply` writes the loops; no virtual calls there.
- The memory report already follows the charges, so it needs nothing per module.
- The per-axis copies, `ex/ey/ez_update_scan`, `update_hx/hy/hz`, `update_psi_e{x,y,z}`, become
  one template on the axis each, as `e_scan<A>` already is.

## Design: materials as definitions loaded into slots

The kernel analogy here is the device tree. Hardware is described as data, a node with a
`compatible` string and properties; a driver declares the strings it matches; at boot each node is
bound to its driver, which probes it and allocates its state. File descriptors are the other half:
a process does not index the kernel's global table, it holds small integers into its own table of
the files it opened.

| kernel | here |
|---|---|
| device-tree node | one material definition: a name, a model, the model's parameters |
| `compatible = "..."` | `model = drude_2cp` |
| driver match table | the list of model types, each with its `name` |
| probe | the model computes its update coefficients for the run's time step |
| per-device state | the slot's polarisation and pole arrays, owned by the slot |
| fd table | the run's slot table; the material map holds slot numbers |
| `module_param` | a parameter table per model: key, member, unit, description |

A definition, one file per material in `materials/`, the catalog:

```ini
; materials/silver_palik_d2cp.ini
[material]
name = silver_palik_d2cp
model = drude_2cp
reference = Vial et al. 2011, fitted to Palik 400 to 1000 nm
eps_inf = 1.4447
plasma_frequency = 1.3280e16
plasma_damping = 9.1269e13

electron_density = 5.541318196761572e28

[critical_point 1]
amplitude = -1.5951
frequency = 8.2749e15
damping = 5.1770e15
phase = 3.1288

[critical_point 2]
amplitude = 0.25261
frequency = 6.1998e15
damping = 5.4126e14
phase = -1.5066
```

These are today's code 101 values. The electron density is read by the hydrodynamic and nonlocal
models; the GNOR diffusion and the nonlocal beta, set from run keys today, become parameters too.

A run names the materials it uses; a `[material]` section in the run file defines or overrides
one by name:

```ini
[materials]
catalog = materials
structure = silver_palik_d2cp
first = air
second = glass

[material glass]
model = dielectric
n = 1.45
```

In the code:

- `materials::Model` is a concept: a `name`, a `Params` struct with its parameter table, a
  `Coeffs` struct, a `State` struct, `probe(params, dt)` and `permittivity(params, omega)`. The
  existing model headers already have `Coeffs`, `compute` and the stretch; they gain `Params` and
  `State` and lose their reads of the global table.
- The parameter table lists each key once, `{"plasma_frequency", &Params::plasma_frequency, ...}`;
  one template parses, validates, prints to the resolved config and documents every model from it.
  The config already does this with `v.field`.
- `Slot` holds a name, `std::variant` of the bound models, `Bound<M>{params, coeffs, state}`, and
  the material's nonlinear tensors. The run's `std::vector<Slot>` has as many entries as materials
  in use, usually two to five.
- The material map holds `std::uint8_t` slot numbers; the stretch dispatch is a `std::visit` on
  the slot, once per stretch, as `with_model` is now.
- chi2 and chi3 become properties of a material, `chi2 = ...` in its definition, so the four
  copies per medium position go.
- The run-dependent entries of today's table, ITO from the carrier density, the graphene stretch,
  the ITO film selected by `middle_layer`, become plain parameters of their definitions; the keys
  that steered them go.

What goes: `codes.h` and its ranges, `ModelParameters`, `coeff_table<M>[500]`, `pc[500]`,
`set_parameters`, the code offsets in geometry and config.

Migration, bit-identical:

1. A script turns every entry `set_parameters` actually sets into a definition file, the doubles
   printed with `%.17g`, named from its comment, the code kept as a comment line. You decide which
   of the historical entries stay in the catalog.
2. Slots replace codes in the update; same coefficients, same order of operations, so the suites
   compare bit for bit once their configs name materials. A script converts the suite configs.
3. The old table and codes are removed.

The one open point: the geometry refers to materials by code in each shape, `material = 101`. With
names, a shape says `material = silver_palik_d2cp` and the renderer asks the slot table for the
number. A shape with no material keeps the structure's.

## Report names

Every report gets comment lines above its header: solver version, precision, polarisation angle,
spectrum range and points, mesh step, time step and steps. The polarisation and the checkpoint
leave the file name where they are metadata; the checkpoint stays in the name of the files written
per checkpoint.

| now | new |
|---|---|
| `Qcoeff_double_<pol>.txt` | `cross_sections.txt` |
| `Qcoeff_NL_double_<pol>.txt` | `cross_sections_harmonic.txt` |
| `Qcoeff_sigma_double_<pol>.txt` | `cross_sections_sigma.txt` |
| `ART-TFSF_double_<pol>_<n>.txt` | `art_tfsf_<n>.txt` |
| `ART-TFSF_NL_double_<pol>_<n>.txt` | `art_tfsf_harmonic_<n>.txt` |
| `ART-TFSF_sigma_double_<pol>.txt` | `art_tfsf_sigma.txt` |
| `ART_<pol>_<n>.txt`, `ART_<pol>_sigma.txt` | `art_<n>.txt`, `art_sigma.txt` |
| `NL_Efficiency_PBC.txt`, `_noPBC.txt`, `_noPBC-ART.txt` | `harmonic_efficiency_periodic.txt`, `harmonic_efficiency.txt`, `harmonic_efficiency_art.txt` |
| `Multipoles_<method>_double_<pol>.txt` | `multipoles_<method>.txt` |
| `Multipoles_<method>_full_current_double_<pol>.txt` | `multipoles_<method>_full_current.txt` |
| `DFT_module_phase.txt` | `source_spectrum.txt` |
| `Signal.txt` | `source_signal.txt` |
| `1D_Ex.txt` ... | `aux_grid_ex.txt` ... |
| `log<rank>.txt` | `run.log`: only root writes it, the number was always 1 |

Found on the way: `DFT_module_phase.txt` has six columns under a four-column header, and
`Signal.txt` has no header; both are fixed with the rename. With `pol` empty the old names became
`Qcoeff_double_.txt`.

## Proposed order

Each step bit-identical in its numbers, checked on all suites, `make regression`, the unit tests
and `check-layers`.

1. Report names and metadata headers. File names change, numbers do not. Done 2026-11-21.
2. Materials as definitions and slots, by the migration above. Done 2026-11-19.
3. Module state into owned structs and the monitor list, module by module, monitors first; the
   wrapper classes and the per-module lines in setup, release and state_file go as each module
   moves.
4. Axis templates for the H update, the CPML faces and the scan wrappers.
5. Namespaces that do not match their directory: `geom`, `h5io`, `kernels`, `fs`, `sm`, the
   one-dimensional simulators and `render_parts`; `inc/app/archive.h` to `inc/io/`.

Decided: the catalog is text files outside the binary; the materials come before the module
state.
