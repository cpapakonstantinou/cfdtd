# Materials

A material is an integer code. The code selects the dispersion model and the parameter set:
`[geometry] material`, `first_medium`, `middle_medium`, `last_medium`, the materials of
`[block.<name>]` sections and the fills of `[anapole]` all take a code.

| Codes | Model | Parameters |
|---|---|---|
| 0 | vacuum | |
| 1 – 99 | dielectric, constant refractive index | `[materials] n<i>` |
| 100 – 149 | Drude + critical points | built in |
| 150 – 199 | Drude | built in |
| 200 – 249 | Lorentz | built in |
| 250 – 299 | modified Lorentz | built in |
| 300 – 349 | hydrodynamic free electrons + critical points | built in, or `[custom_material]` |
| 350 – 399 | Drude + one Lorentz pole, for the two-temperature and hydrodynamic electron gas | built in, or `[custom_material]` |
| 400 – 449 | magnetic Lorentz | built in |

## Dielectrics

```ini
[geometry]
middle_medium = 4

[materials]
n4 = 1.45
```

`n<i>` is required for every dielectric code the run uses.

## Built-in dispersive materials

| Code | Material | Model |
|---|---|---|
| 100 | gold | Drude + 2 critical points |
| 101 | silver | Drude + 2 critical points |
| 102 | graphene (300 – 800 nm) | Drude + 4 critical points; `[materials] stretch_x` scales it |
| 103 | gold-silver alloy AuAg 0.5 (270 – 1200 nm) | Drude + critical points |
| 104, 105 | gold, fits to Johnson and Christy (188 – 1937 nm) | Drude + critical points |
| 106 | chromium, Sytchkova (262 – 2500 nm) | Drude + critical points |
| 150 | gold, McPeak (1 – 1.7 µm) | Drude |
| 151 – 156 | indium tin oxide (ITO) fits; the plasma frequency follows `[materials] carrier_density` (151, 153, 155) or `carrier_density_back` (152, 154, 156) | Drude |
| 200 | alumina (0.2 – 1 µm) | Lorentz |
| 201 | fused silica (0.2 – 1.5 µm) | Lorentz |
| 250, 252 | silicon | modified Lorentz |
| 300 | gold | hydrodynamic |
| 301 | silver | hydrodynamic |
| 302 | the material of `[custom_material]` | hydrodynamic |
| 304 | gold, Johnson and Christy | hydrodynamic |
| 350 | thin-film ITO; the parameter set follows `[geometry] middle_thickness` ([ITO](ITOmodel.md)) | Drude + Lorentz |
| 351 | the material of `[custom_material]` | Drude + Lorentz |

The parameter sets are in `src/materials/parameters.cc`. Several slots that are not in this table
(202, 203, 251, 400 and the other ITO variants) hold fits whose labels in the source are historical;
check the values there before using one.

## A custom material: `[custom_material]`

Codes 302 and 351 take their parameters from the file:

```ini
[custom_material]
eps_inf = 1.0
plasma_frequency = 1.37e16     ; [rad/s]
plasma_damping = 1.0e14        ; [rad/s]
cp_frequency = 4.0e15          ; critical point 1 (the Lorentz pole of code 351)
cp_damping = 8.0e14
cp_amplitude = 1.0
cp_phase = -0.785
cp_frequency_2 = 5.6e15        ; critical point 2 (code 302 only)
cp_damping_2 = 2.0e15
cp_amplitude_2 = 1.2
cp_phase_2 = -0.785
fermi_velocity = 1.0           ; code 302: scales the Fermi velocity
electron_density = 5.9e28      ; code 302: equilibrium density [1/m^3]
```

The numbers above only show the format.

## Models that act on a material

| Model | Selected by | Needs |
|---|---|---|
| nonlocal, GNOR | `[electrons] model = nonlocal / gnor` | a Drude + critical points material (100 – 149); [Nonlocal Plasmonics](NonlocalPlasmonics.md) |
| hydrodynamic | `[electrons] model = hydrodynamic` | a hydrodynamic material (300 – 349); [Hydrodynamics](Hydrodynamics.md) |
| two-temperature, pump-probe | `[temperature] model` | a Drude + Lorentz material (350 – 399); [ITO](ITOmodel.md) |
| χ² and χ³ | `[nonlinear] shg / thg` | any material; [Nonlinear Tensors](NonlinearTensors.md) |
