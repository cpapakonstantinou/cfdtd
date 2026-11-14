# Where the shared state goes

`fdtd::state` was the old global state; it was removed on 2026-11-13. This page records where each of its variables went (names in the tables are the names of that time; most were renamed to snake_case since). Variables used by one module live in that module (`inc/<module>/state.h`). This page covers the rest: values that several modules read. Each one gets one owner, chosen by **when it is written**.

| Layer | Written | Owner | Read as |
|---|---|---|---|
| Run configuration | once, in the configure phases | `app::RunConfig` (`inc/app/run_config.h`), sealed by the driver after rank initialisation | `run_config().grid.dx` (const) |
| Rank decomposition | once per rank, by the decomposition | `grid::Subdomain` | `subdomain.…` |
| Mesh products | per sweep, when the mesh and materials are built | the module that builds or updates them (`materials`, `nonlinear`, `hydro`, `thermal`, `monitors`) | passed in, or the module's object |
| Run state | while the run is going | `app::` driver context (step, sweep, per-sweep polarisation) | `StepTime` / `SampleTime` arguments |

Why split it this way:
- A value that is fixed for the whole run is configuration. Making it const after configure turns "a global anyone can change" into a contract.
- The configuration is resolved once, so flags such as `second_harmonic` or `hydrodynamics` are tested where the driver assembles the run (the source and monitor lists already work like this), not inside the kernels.
- Values that are per rank, per mesh or per step are not configuration. Putting them in `RunConfig` would recreate `fdtd::state` under a new name.

## 1. Run configuration: done (2026-11-07)

154 variables, in nine groups: `time`, `grid`, `source`, `materials`, `nonlinear`, `thermal`, `geometry`, `analysis`, `output`.

**Access:**
- `.cc` files cache `const RunConfig &run_cfg = run_config();` once, at file level.
- Headers call `fdtd::app::run_config()`, which is inline.
- Only the configure phases take write access, through `edit_run_config()`. These are `apply_sim_config`, `app::setup_domain`, `app::initialise_run`, `Driver::load_configuration`, `Driver::initialize_ranks` and `assign_output_paths`. Any other write fails to compile, and a late call to `edit_run_config()` throws.

**Next:**
- [x] `pids` (`INI_DATA_STRUCT`) retired (2026-11-07): `apply_sim_config` fills `RunConfig` directly. Its 90 fields joined the groups; where a field duplicated a `RunConfig` value (`WIDTH`, `pbc_x`, `save_output`, ...) the two are one member. The per-rank cell counts of the decomposition went to `grid::cells_per_rank_x/y/z`.
- [x] The 46 shape flags (34 in `RunConfig`, 12 in `fdtd::geometry`), each true for exactly one `antenna_type`, are gone (2026-11-07): code tests `run_cfg.geometry.is(shape::sphere)`. The two derived flags stay members (`grating`, `meta_butterfly`: a condition over several shapes and settings).
- [x] Core modules no longer read `RunConfig` at all (2026-11-13): the application hands each one a typed input ([Architecture](Architecture.md)).
- [x] Upper-case legacy names renamed (2026-11-13). A few long-hand ones remain (`seperate_linear_nonlinear`, `nr_cells_no_cmpl_*`).

## 2. Rank decomposition (69 variables): to `grid::Subdomain`

| Variables | Meaning |
|---|---|
| `startindices`, `localsize`, `rank_this`, `global_procs_size` | rank's place in the global grid |
| `nea_/far_cpml_x/y/z` | rank touches a CPML face |
| `IPW1`, `JPW1`, `KPW1`, `ipw_2` … `kpw_max`, `j_for_uni`, `has_plane_wave` | local TF/SF box |
| `i/j/k_qsca_*`, `i/j/k_qtot_*` | local flux boxes |
| `j_art_1/2` | local ART planes |
| `animation_*`, `i/j/k_animation*` | local animation planes |
| the halo MPI datatypes, the output communicators | rank's communication handles |

Several of these duplicated fields that `Subdomain` already has (`subdomain.ipw_1` was copied into `IPW1`). The globals go once the readers use the `Subdomain` field.

- [x] Done (2026-11-07), 54 variables: the TF/SF box indices (`IPW1` ... `kpw_max`), the ART plane indices, the flux-box indices (`[ijk]_q{sca,tot}_{min,max,nea,far}`), the animation flags and indices (`animation_*` are `subdomain.plots_*()`, the indices `subdomain.[ijk]_animation*`) and the CPML-face flags (`subdomain.nea_cpml_x()` ..., from `which_type_of_process_*`). The copies were conditional (left at 0 when the feature was off); an instrumented run compared every copy with its field on the regression layouts, and the differences were only on features that were off, whose readers do not run then.
- [x] `startindices`, `localsize`, `j_for_uni`, `has_plane_wave` are `Subdomain` members (`output_start`, `output_size`, `j_for_uni`, `has_plane_wave`), set where they were set (2026-11-07).
- [x] `rank_this` / `global_procs_size` are `grid::world_rank` / `world_size` (set in `main`; `subdomain.rank` exists only after the decomposition); the halo datatypes and the output communicators are `grid::` state (2026-11-07). Section 2 is complete.

## 3. Mesh products (138 variables): to the module that builds them

Done (2026-11-07) for the owners that need no redesign: the temperatures to `thermal`; the chi tables, nonlinear flags and cell maps to `nonlinear` (new `inc/materials/nonlinear/state.h`); the electron-gas arrays and constants to `hydro`. 61 variables, ownership only (uses qualified from the build errors). Left: the material coefficients and ADE state (they become `Coeffs` / `State` of the model types), the material maps, the signal normalisation, `energy_density`, the movie files.

| Variables | Owner |
|---|---|
| `matind`, `matind_x/y/z[_full]` | `materials::MaterialMap`, built by `geometry` |
| `ce_1`, `ce_2`, `ce_*_d`, `ce_*_l`, `epre`, `sig_m`, `n_refr*`, the Drude/2CP/Lorentz parameters and coefficients (`b_cp_*`, `a_cp_*`, `c_cp_*`, `aj`, `bj`, `alpha_j_*`, `beta_j`, `nonlin_drude_*`, …) | `materials`: one `Coeffs` per model (§3 step 4 of the roadmap) |
| `Px/Py/Pz`, `p*_prev`, `p*_cp[_prev]`, `jx` … `jz_nl` | `materials`: the model `State` (to be packed per stretch) |
| `CHI2`, `chi_2_*`, `chi_3_*`, `chi*_dl`, `iso_chi_*`, `nl`, `nl_dl` | `nonlinear` |
| `ex_2` … `ez_2`, `epsinf_*`, `grad_div_*`, `p*_nl[_2]`, `p*_cp_2[_nl]`, `non_loc`, `beta_nl`, `diff_coeff`, `gnor_*`, electron-gas constants | `hydro` |
| `electron_temp[_tf]`, `lattice_temp`, `temp_save` | `thermal` |
| `magnetization_x/y/z` | `materials` (magnetic Lorentz) |
| `OMEGA`, `amplitude_normalization`, `dft_normalization` | the signal (`sources`) and the analysis frequencies |
| `energy_density` | `monitors` |
| `ex_file`, `ez_file`, `cstring` | driver (movie output) or removed |

## 4. Run state (35 variables): to the driver, or passed as arguments

| Variables | Owner |
|---|---|
| `t`, `t_1`, `timestep`, `prgstart` | `StepTime` / driver context |
| `sweep`, `xpol`, `zpol`, `pol`, `nr_max_time_steps`, `input_power`, `tt` | per-sweep state in the driver context; the signal object for `tt` and `input_power` |
| `chx`, `chy`, `chz` | `update` (H coefficients) |
| `focus_*`, `x/y/z_tfsf_n/f`, `offset_delay` | the superposition source (roadmap §3 step 5) |
| `is_chi_*_a_*`, `plasma_freq_add_term`, `nonlocal_barrier_int`, `use_backwards_differences` | `nonlinear` / `materials` / `hydro` |
| `dfg_optim`, `count_art`, `geometric_cross_section`, `width`, `failure_type` | driver context, ART monitor, geometry |

## Status (2026-11-13)

`fdtd::state` is gone, and so is `inc/app/globals.h` (the header that only gathered the module headers): every file includes the headers it uses. Each source file names what it reads with using-declarations (`using fdtd::fields::Ex;`, `using fdtd::materials::ce_2;`), so the dependencies of a file are listed at its top. The objects that replaced groups of globals: `fields::Fields`, `grid::Subdomain`, `materials::MaterialMap` (the material maps), `sources::Signal` (the time signal, held by every source), `WaveFrame` (the anchor of the superposed waves, held by the superposition source and the 1D simulators). No module reads the run's time or step: the time loop passes `StepTime` / `SampleTime` or the step. `inc/app/globals.h` still gathers the module headers (include hygiene is open).

## Status (2026-11-11, before the removal)

All four layers are assigned: `inc/app/globals.h` declares no variables. The run state is `app::RunState` (`inc/app/run_state.h`: `t`, `t_1`, `timestep`, `sweep`, `nr_max_time_steps`, `prgstart`, `failure_type`, `dfg_optim`, the trace files). The mesh products left in §3 (material tables, maps, P, poles, currents) are `fdtd::materials` state; the per-sweep polarisation, signal normalisation and superposition geometry of §4 are `fdtd::sources` state. Still to do: pass the time as `StepTime` / `SampleTime` instead of reading `run_state()`, and turn the module state files into objects.
