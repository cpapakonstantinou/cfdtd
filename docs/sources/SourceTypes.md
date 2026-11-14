# Source Types

A source is a spatial profile and a time signal:

```ini
[source]
type = plane_wave        ; plane_wave | wave_superposition | gaussian_beam | dipole | mode
signal = gaussian        ; the time signal; its keys are in [signal.gaussian]
peak_field = 1           ; [V/m]; matters for nonlinear runs
```

The keys of a type are in the section named after it, and the keys of a signal in
`[signal.<name>]`. Setting a key of a type or signal the run does not select is an error.

## Spatial profiles

| `type` | Injects | Keys |
|---|---|---|
| `plane_wave` (default) | a plane wave along +y on the TF/SF box | `[plane_wave]` |
| `wave_superposition` | any number of plane waves at arbitrary angles | `[wave.<name>]`, `[wave_superposition]`; [Plane Wave Superposition](PlaneWaveSuperposition.md) |
| `gaussian_beam` | a focused beam, built from plane waves | `[gaussian_beam]`; [Gaussian Beam](GaussianBeam.md) |
| `dipole` | point dipoles | `[dipole.<name>]` |
| `mode` | a waveguide mode computed with MPB | `[mode]`; [Mode Source](ModeSource.md) |

### Plane wave

```ini
[plane_wave]
polarization_axis = z      ; x | z | rhcp | lhcp
polarization_angle = 0     ; [deg] turns E from polarization_axis toward the other transverse axis
```

The wave travels along +y from the first medium. It is injected on the faces of the TF/SF box
(`[boxes] tfsf_*`); a 1D auxiliary grid carries the incident field through the layers of
`[geometry]`, so the layers and the grid dispersion are accounted for. With periodic x and z the wave
is injected on one plane. `rhcp` and `lhcp` give circular polarisation.

Oblique incidence is `type = wave_superposition` with one wave.

### Dipoles

One section per dipole; positions are in cells, the amplitudes multiply the signal:

```ini
[source]
type = dipole

[dipole.a]
x = 100
y = 100
z = 100
e_x = 0
e_y = 0
e_z = 1
```

## Time signals

| `signal` | Shape | Keys of `[signal.<name>]` |
|---|---|---|
| `raised_cosine` (default) | a short broadband pulse that starts and ends at exactly zero; use it for linear spectra | `max_frequency` [Hz] |
| `raised_cosine_derivative` | the time derivative of the raised cosine (no DC content) | `max_frequency` |
| `rectangular` | a top-hat pulse in time, as long as the raised cosine of the same `max_frequency`; not well tested | `max_frequency` |
| `sinc` | a windowed sinc on a 520 nm carrier, flat over about 350 – 1000 nm (Hassan and Calà Lesina, Opt. Express 30, 19557, 2022); the band is fixed in the code | none |
| `gaussian` | a Gaussian envelope on a carrier; use it for nonlinear runs and for oblique incidence on periodic structures | `center_wavelength`, `width` (FWHM of the field) [s], `chirp`, `chirp_quadratic` |
| `gaussian_alt` | two Gaussian pulses on two carriers; cannot be chirped | `center_wavelength`, `width` |
| `cw_pulse` | a carrier at the centre of `[spectrum]`, switched on and off smoothly | `periods` |
| `continuous_wave` | a carrier at the centre of `[spectrum]` that stays on | `settle_periods`: the DFT starts after this many periods |
| `constant` | a field that ramps up and holds | none |

The run stops before it starts if the time steps do not cover the signal;
[Signal and Time Steps](../tutorial/SignalSettings.md) shows how to choose them with `fdtd signal`.
`reports/Signal.txt` holds the signal the run used.
