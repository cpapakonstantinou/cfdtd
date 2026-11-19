# Migrating Old Input Files

Before the INI format a run was described by `pphinfoini.json` or `pphinfoini.txt` plus side files
in `INIDEF/` (`pphmatini.txt`, `plane_waves.json`, `chi2/`, ...). Those files are no longer read.
This page maps their keys and codes to the INI keys.

`tools/migrate_config_keys.py run.ini` rewrites an INI file written for an earlier state of the INI
schema (renamed keys, the generic `[geometry]` dimensions, rotation conventions) in place and
reports what it changed or dropped. `tools/migrate_materials.py run.ini` then turns the material
codes of such a file into material names, [Materials](../models/Materials.md).

Two conventions changed with the format:

- `[grid] nx`, `ny`, `nz` count interior cells; the old `Domain Size` included the PML on both sides
  (`nx = Domain Size - 2 * PML Box`).
- The spectrum keys mean the same in every simulation. The solver no longer replaces the first two
  points by the pump and its harmonic for hydrodynamic or temperature runs; ask for that pair with
  `lambda_min = lambda_max = <pump>`, `points = 1`, `harmonic_band = second`.

## Keys

| old (pphinfoini.json / .txt / side files) | new INI key |
|---|---|
| Domain Size [nx,ny,nz] | `[grid] nx`, `ny`, `nz` (required) |
| Space Step / Step Size [dx,dy,dz] | `[grid] dx`, `dy`, `dz` (required) |
| PML Box | `[grid] pml_x`, `pml_y`, `pml_z` |
| Number of Time Steps | `[grid] time_steps` |
| Center Position | `[grid] origin_x`, `origin_y`, `origin_z` |
| Periodic X / Z | `[grid] periodic_x`, `periodic_z` |
| Staircasing Type | `[grid] staircase = true/false` |
| cfl_coeff (hardcoded 0.5) | `[grid] courant` |
| Domain Decomposition | `[parallel] ranks_x`, `ranks_y`, `ranks_z` |
| Chunks StartIndex X/Y/Z | `[parallel] chunks_x`, `chunks_y`, `chunks_z` |
| TFSF Box / Scattering Box / Total Box | `[boxes] tfsf_*`, `scattering_*`, `absorption_*` |
| Source Type / TFSF 1/2/3/4 | `[source] type = plane_wave / dipole / gaussian_beam / wave_superposition / mode`; the keys of a type are in the section of the same name |
| Signal Type 1..11 | `[source] signal = raised_cosine / raised_cosine_derivative / continuous_wave / cw_pulse / rectangular / gaussian / gaussian_alt / constant / mode / sinc`; the keys of a signal are in `[signal.<name>]` |
| Center Wavelength, Pulse Width, Maximum Frequency | `[signal.gaussian]` (also `gaussian_alt`) `center_wavelength`, `width`, `chirp`, `chirp_quadratic`; `[signal.raised_cosine]` (also `raised_cosine_derivative`, `rectangular`) `max_frequency` |
| Polarization Angle (2) | `[plane_wave] polarization_axis`, `polarization_angle` (one polarisation; a second one is a second run) |
| Maximum Field | `[source] peak_field` |
| Polarization, circ_LH / circ_RH / circ_time_shift | `[plane_wave] polarization_axis = x / z / rhcp / lhcp`, `polarization_angle` (degrees from that axis toward the other transverse axis); `propagation_axis = y`, `propagation_angle = 0` |
| num_periods_cw, start_cw, PW_factor (hardcoded) | `[signal.cw_pulse] periods`, `[signal.continuous_wave] settle_periods`, `[source] aux_length_factor` |
| 1DSim Number of Time Steps, Enforce 1D Simulation | `[source] aux_steps`, `[wave_superposition] layered_simulation` |
| gaussian_beam.txt | `[gaussian_beam] theta_steps`, `phi_steps`, `theta_max`, `phi_max`, `waist_factor`, `focal_factor` |
| Number/Minimum/Maximum Wavelength (meaning changed with the run type!) | `[spectrum] lambda_min`, `lambda_max` (required), `points`, `spacing`, `harmonic_band` |
| DFT Plot 1/2/3/4/5, Full Domain DFT | `[output] dft_yz`, `dft_xz`, `dft_xy` (any combination; 4 = all three), `dft_magnetic` (5 = all three + H), `volume_dft` |
| Movie XY/XZ/YZ/3D | `[output] frames_xy`, `frames_xz`, `frames_yz`, `frames_3d` |
| DFT Phase Normalization, Sigma Absorption | `[output] phase_normalize`, `sigma_absorption` |
| (new) | `[output] xdmf` (default true): a `.xmf` next to every HDF5 output for ParaView / VisIt; complex fields are stored as `<name>_re` / `<name>_im` |
| Exact Multipoles | `[multipoles] exact`; new: `lwa` (long-wavelength moments, any material model) |
| Multipoles [r, ntheta, nphi] | `[multipoles] sphere_radius`, `sphere_theta`, `sphere_phi` |
| Antenna Type 0..1001 / Antenna Medium | `[geometry] shape = none / sphere / dipole / bowtie / ...` (names below), `material`, a material name |
| First / Middle / Last Medium, Middle Thickness | `[geometry] first_medium`, `middle_medium`, `last_medium`, material names; `middle_thickness` |
| Antenna Shift | `[geometry] shift_y` |
| WIDTH, LENGTH, THICKNESS, LENGTH2, GAP, LENGTH1_L, LENGTH2_B | keys of the shape's own section, named for the shape ([Shapes](../geometry/Shapes.md)); `Height` is removed (read by nothing) |
| ald_layer, afm_size_x/z (hardcoded) | removed (read by nothing) |
| pphmatini.txt NREFR1..99 | a `[material.<name>]` of model `dielectric` with `n` per medium |
| pphmatini.txt carrier/gamma/stretch/sheets | `carrier_density` of the ITO materials, `sheet_scale` of graphene (`stretch_z`, `sheet_spacing`, and the LIDAR-only `damping_*` / `extra_layer_*` are removed) |
| pphREADINMATERIAL.txt | a `[material.<name>]` of model `hydrodynamic` or `drude_lorentz` with its poles |
| pphlidar.txt | removed with the `lidar_*` shapes (2026-11-13) |
| DFG_OPTIM 0..6 | `[nonlinear] shg = true / false`, `thg = true / false`; the materials with a tensor are nonlinear (the DFG/SFG two-sweep modes are removed: run each pump separately and mix in post-processing) |
| pphnlsini.txt (4 tensor names + approach) | `chi2 =` / `chi3 =` of the material definitions; `[nonlinear] update = field / polarization` |
| INIDEF/chi2/*.txt, INIDEF/chi3/*.txt | `[chi2.<name>]` / `[chi3.<name>]` sections, in the run file or in the catalog, `materials/chi2/`, `materials/chi3/` ([Nonlinear Tensors](../models/NonlinearTensors.md)) |
| DFG_THZ_FREQ, DFG_DEEP_STEP, GRAPHENE_OPTIM | `rotation` [deg] and `corner_radius` [m] of the shape's section (what the shapes read the two DFG values as); `[nonlinear] graphene_optimization` |
| Advanced Model 0..2023 | `[electrons] model`, `[temperature] model`, `[hydrodynamic]` (table below) |
| Diffusion Coefficient, Fermi velocity | `diffusion` of the material; `[electrons] fermi_velocity`, a factor on the materials' `fermi_velocity` |
| NumRK.txt, PumpProbeDelay.txt | `[two_temperature] rk_stages`, `[pump_probe] delay` |
| Probes.txt | one `[probe.<name>]` section with `x`, `y`, `z` per probe |
| ./INIDEF, ./OUTPUT (hardcoded) | `[paths] input_dir` (data files: mode profiles, electron temperature), `output_dir`; CLI `-i`, `-o` |
| monitors.json | one `[monitor.<name>]` section per box: `type = frequency / time`, `center_x/y/z`, `size_x/y/z` [cells], `step` |
| Dipole_Source_Locs.txt, Dipole_Amps.txt | one `[dipole.<name>]` section per dipole: `x`, `y`, `z` [cells], `e_x/y/z` |
| plane_waves.json | one `[wave.<name>]` section per wave: `direction_x/y/z`, `e_x/y/z`, `h_x/y/z`; `[wave_superposition] derivative`, `dispersion_correction`, `layered_simulation` |
| dfts.json, geometry.json, adjoint gradients | removed (no JSON input; use `[monitor.<name>] type = frequency`) |
| 23.5 h wall-clock stop (hardcoded) | `[run] wallclock_limit_hours` (0 = never) |

## "Advanced Model" codes

| old code | new keys |
|---|---|
| 0 | `[electrons] model = local` (default) |
| 1 / 2 | `[electrons] model = nonlocal` / `gnor` |
| 3 | `[electrons] model = hydrodynamic` (`[hydrodynamic] pressure`, `convection`, `magnetic` default on) |
| 4 .. 9 | hydrodynamic with `[hydrodynamic] pressure` / `convection` / `magnetic` switched off as needed (4: no pressure; 5: magnetic only; 6: none; 7: convection only; 8: no magnetic; 9: no convection) |
| 10 / 11 / 12 | like 6 / pressure only / 3 |
| 13 / 14 / 15 | like 3: the density update, the SHG source terms and the critical-point switch of these codes are not wired in the solver (roadmap §7), so their keys are gone |
| 200 / 201 | 4 / 3 + `[hydrodynamic] split_nonlinear = true` |
| 202 / 2020 | `[temperature] model = two_temperature` (2020: `[two_temperature] save_frames = true`) |
| 2021 | `[temperature] model = load` (reads `<input_dir>/electron_temperature.h5`) |
| 2022 / 2023 | `[temperature] model = pump` / `probe` (`[pump_probe] delay`) |
| 203 | hydrodynamic with `pressure = convection = magnetic = false`, `[temperature] model = two_temperature`, `[two_temperature] save_frames = true` |

`[electrons] model = hydrodynamic` requires a hydrodynamic material (codes 300 to 349).

## "Antenna Type" codes

Shapes that are not in this table were removed.

| code | name | code | name | code | name |
|---|---|---|---|---|---|
| 0 | `none` | 1 | `sphere` | 5 | `dipole` |
| 6 | `bowtie` | 12 | `cylinder` | 13 | `two_spheres` |
| 17 | `l_antenna` | 19 | `hole` | 50 | `lateral_dipole` |
| 54 | `round_butterfly` | 55 | `round_dipole` | 56 | `one_side_round_dipole` |
| 57 | `round_tripole` | 58 | `round_hybrid` | 59 | `x_shape` |
| 67 | `sphere_dimer` | 68 | `sphere_trimer` | 69 | `sphere_quadrumer` |
| 86 | `ellipse` | 100 | `sphere_shell` | 103 | `triangle` |
| 104 | `rounded_triangle` | 105 | `rounded_monopole` | 106 | `u_shape` |
| 108 | `cylinder_3d` | 109 | `h_shape` | 110 | `v_shape` |
| 111 | `l_shape` | 112 | `two_films_spacer` | 113 | `rounded_dipole` |
| 114 | `rounded_dipole_offset_planes` | 115 | `cylinder_shell` | 2001 | `anapole` |
