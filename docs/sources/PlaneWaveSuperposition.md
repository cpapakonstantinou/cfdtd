# Plane Wave Superposition

## Concept

`[source] type = wave_superposition` injects a list of plane waves, each with its own direction and
field amplitudes, and superimposes them on the TF/SF box. One wave gives oblique incidence; many
waves build structured light. `[source] type = gaussian_beam` builds such a list for a focused beam
([Gaussian Beam](GaussianBeam.md)).

## How to use

Each wave is one `[wave.<name>]` section:

| Key | Meaning |
|---|---|
| `direction_x`, `direction_y`, `direction_z` | unit vector of propagation |
| `e_x`, `e_y`, `e_z` | electric field amplitude; multiplies the signal |
| `h_x`, `h_y`, `h_z` | magnetic field amplitude: `direction × e / η`, with η = 376.73 Ω in vacuum |

and three switches apply to all waves:

| `[wave_superposition]` key | Default | Meaning |
|---|---|---|
| `derivative` | false | inject the time derivative of the signal (the focused-beam construction needs it) |
| `dispersion_correction` | true | correct the phase velocity of each wave for grid dispersion ([Numerical Dispersion Correction](NumericalDispersionCorrection.md)) |
| `layered_simulation` | false | run the 1D simulation of each wave even when the three layers are one medium ([Plane Wave Superposition with Media](Oblique1DSimulations.md)) |

When the layers of `[geometry]` differ, the 1D simulations run automatically.

## Example

Two waves at ±45° in the x-y plane, polarised along z:

```ini
[grid]
nx = 200
ny = 320
nz = 200
dx = 50e-9
dy = 50e-9
dz = 50e-9
time_steps = 1500

[parallel]
ranks_x = 2
ranks_y = 2
ranks_z = 2

[boxes]
tfsf_x = 20
tfsf_y = 20
tfsf_z = 20

[spectrum]
lambda_min = 0.5e-6
lambda_max = 1.5e-6
points = 31

[source]
type = wave_superposition
signal = gaussian

[signal.gaussian]
center_wavelength = 1e-6
width = 15e-15

[wave.left]                 ; travelling toward -x and +y, 45 deg
direction_x = -0.7071067811865476
direction_y = 0.7071067811865476
direction_z = 0
e_x = 0
e_y = 0
e_z = 1
h_x = 0.00265442
h_y = 0.00265442
h_z = 0

[wave.right]                ; travelling toward +x and +y, 45 deg
direction_x = 0.7071067811865476
direction_y = 0.7071067811865476
direction_z = 0
e_x = 0
e_y = 0
e_z = 1
h_x = 0.00265442
h_y = -0.00265442
h_z = 0

[output]
frames_xy = 15
frames_xz = 15
```

## Notes

- **`derivative`:** the plane-wave expansion of a focused beam contains the time derivative of the
  signal; set the key for a hand-built superposition of that kind.
- **`dispersion_correction`:** the correction uses the plane-wave dispersion relation at the carrier
  of a Gaussian signal, so it is exact at the centre frequency only and still imperfect on very
  coarse grids.
- **Restart:** a run with layered 1D simulations cannot be continued from a state file
  ([Stopping and Restarting](../tutorial/Restart.md)).
