# Materials

A material is a definition: a name, an update model and the model's parameters, written as INI
sections. The run names the materials it uses; `[geometry] material`, `first_medium`,
`middle_medium`, `last_medium`, `shell_material`, the `material` of each `[block.<name>]` and the
fills of `[anapole]` all take a name. A definition is found in the run file first, then in the
catalog, a directory of files that the solver reads at start-up. Adding or changing a material is
editing a text file; nothing is compiled.

```ini
[geometry]
material = silver_d2cp      ; catalog file materials/silver_d2cp.ini
last_medium = glass         ; defined below

[material.glass]
model = dielectric
n = 1.45
```

`vacuum` is built in. It is the material of every cell nothing is drawn in, and the default of every
key above.

## The catalog

The solver looks for `[material.<name>]` in `<name>.ini` of each catalog directory, in this order:

1. the directories of `[materials] catalog`, separated by `:`, relative to the run file;
2. the directories of the environment variable `FDTD_MATERIALS`;
3. `materials/` next to the executable, the catalog shipped with the solver.

The first file found is used. The tensors named by `chi2 =` and `chi3 =` are found the same way, in
`chi2/<name>.ini` and `chi3/<name>.ini`. The resolved configuration of a run,
`resolved_config.ini`, holds every material and tensor it used, so the run can be repeated without
the catalog.

## A definition

```ini
[material.silver_d2cp]
model = drude_cp
reference = Drude + 2 critical points fit of silver, Vial et al. 2011, to Palik 400 to 1000 nm
eps_inf = 1.4447
plasma_frequency = 1.328e16     ; [rad/s]
plasma_damping = 9.1269e13      ; [rad/s]

[material.silver_d2cp.pole.1]
amplitude = -1.5951
frequency = 8.2749e15           ; [rad/s]
damping = 5.177e15              ; [rad/s]
phase = 3.1288                  ; [rad]

[material.silver_d2cp.pole.2]
amplitude = 0.25261
frequency = 6.1998e15
damping = 5.4126e14
phase = -1.5066
```

Names are lower case, with digits, `_` and `-`. A key the model does not take is an error, as is
a model without a required key. `reference` is for the reader; the solver does not use it.

`extends = <name>` starts from another material and changes what it lists, key by key and pole by
pole:

```ini
[material.silver_rough]
extends = silver_d2cp
plasma_damping = 1.5e14
```

Any key can be set from the command line as well, `-s material.silver_d2cp.plasma_damping=1.5e14`,
or a material defined there entirely, `-s material.d15.model=dielectric -s material.d15.n=1.5`.

## Models and their keys

| model | keys | poles |
|---|---|---|
| `dielectric` | `n` or `eps`, `conductivity` [S/m] | none |
| `drude` | `eps_inf`, `plasma_frequency` or `carrier_density` [m⁻³] with `effective_mass` [electron masses], `plasma_damping` | none |
| `drude_cp` | `eps_inf`, `plasma_frequency`, `plasma_damping`; `sheet_thickness` [m] and `sheet_scale` for a sheet drawn one cell thick | 1 to 4: `amplitude`, `frequency`, `damping`, `phase` |
| `lorentz` | `eps_inf`, `eps_static`, `resonance_frequency`, `damping` | none |
| `modified_lorentz` | `eps_inf` | 1 to 4: `amplitude`, `frequency`, `damping`, `numerator` [rad/s] |
| `hydrodynamic` | as `drude_cp`, and `effective_mass`, 0 = from the density and the plasma frequency | as `drude_cp` |
| `drude_lorentz` | `eps_inf`, `plasma_frequency`, `plasma_damping` | 1: `amplitude`, `frequency`, `damping` |
| `magnetic_lorentz` | `eps_inf`, `eps_static`, `resonance_frequency`, `damping`, `magnetic_coupling` | none |

Frequencies and dampings are angular, in rad/s. Every model takes `conductivity`. The free-electron
models, `drude`, `drude_cp`, `hydrodynamic` and `drude_lorentz`, also take the electron gas that the
nonlocal and hydrodynamic models read: `electron_density` [m⁻³], `fermi_velocity` [m/s] and
`diffusion` [m²/s], the GNOR diffusion. `[electrons] fermi_velocity` is a factor on every material's
`fermi_velocity`.

## Nonlinear materials

A material is nonlinear when its definition names a tensor:

```ini
[material.gold_kerr]
extends = gold_d2cp
chi3 = au              ; catalog chi3/au.ini, or a [chi3.au] section in the run file

[chi2.my_crystal]      ; d11 ... d36, d_mn in m/V; the solver uses chi2 = 2 d
d22 = 2.1e-12
d31 = 0.15e-12
```

`[nonlinear] shg = true` loads the `chi2` tensors, `thg = true` the `chi3` tensors, and every cell of
a material with a non-zero tensor takes part, wherever the material is drawn. A structure and a
substrate of the same metal where only the structure should radiate use two materials, one
extending the other with the tensor. [Nonlinear Tensors](NonlinearTensors.md) lists the components
and the catalog tensors.

## The catalog shipped with the solver

The files in `materials/` are the parameter sets of the former built-in table, with the values it
held. Each file names the former material code in a comment.

| name | model | former code |
|---|---|---|
| `gold_d2cp`, `silver_d2cp` | Drude + 2 critical points | 100, 101 |
| `graphene_d4cp` | Drude + 4 critical points, a 0.34 nm sheet | 102 |
| `gold_silver_alloy_d2cp` | Drude + 2 critical points, 270 to 1200 nm | 103 |
| `gold_jc_d2cp`, `gold_jc_d4cp` | fits to Johnson and Christy, 188 to 1937 nm | 104, 105 |
| `chromium_d2cp` | Sytchkova, 262 to 2500 nm | 106 |
| `gold_mcpeak_drude` | Drude, McPeak, 1 to 1.7 µm | 150 |
| `ito_lidar_drude`, `ito_10p5_drude`, `ito_2p6_drude` | Drude from the carrier density | 151 to 156 |
| `alumina_lorentz`, `silica_lorentz`, `gaas_lorentz`, `lorentz_test_600thz` | Lorentz | 200 to 203 |
| `silicon_modified_lorentz`, `gaas_modified_lorentz` | modified Lorentz | 250, 251 |
| `gold_hydrodynamic`, `silver_hydrodynamic`, `gold_jc_hydrodynamic` | hydrodynamic | 300, 301, 304 |
| `ito_film_thin`, `ito_film_mid`, `ito_film_thick` | Drude-Lorentz, films under 50 nm, 50 to 300 nm, over 300 nm | 350 |
| `magnetic_lorentz_test` | magnetic Lorentz test medium | 400 |

`tools/migrate_materials.py` rewrites a configuration with material codes to these names.

## Models that act on a material

| Model | Selected by | Needs |
|---|---|---|
| nonlocal, GNOR | `[electrons] model = nonlocal / gnor` | a `drude_cp` material with `fermi_velocity`, `diffusion`; [Nonlocal Plasmonics](NonlocalPlasmonics.md) |
| hydrodynamic | `[electrons] model = hydrodynamic` | a `hydrodynamic` material; [Hydrodynamics](Hydrodynamics.md) |
| two-temperature, pump-probe | `[temperature] model` | a `drude_lorentz` material, the first one is heated; [ITO](ITOmodel.md) |
| χ² and χ³ | `[nonlinear] shg / thg` | a material with a tensor, above |
