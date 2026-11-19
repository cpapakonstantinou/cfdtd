# The Configuration File

A run is described by one INI file. `fdtd template` prints every key with its default and a
one-line description; this page explains how the sections fit together.

## Rules

- `[section]` headers, `key = value` lines, `;` or `#` starts a comment.
- Lengths are in metres, times in seconds, angles in degrees, box and probe positions in cells.
- Every key is optional unless the template marks it *(required)*: `[grid] nx, ny, nz, dx, dy, dz`
  and `[spectrum] lambda_min, lambda_max`.
- An unknown or misspelled key is an error, and so is a key in the section of an option the run
  does not select (a `[sphere]` key with `shape = dipole`, a `[hydrodynamic]` key with
  `[electrons] model = local`). The message names the key.
- `-s section.key=value` on the command line overrides a key; `-o` and `-i` override
  `[paths] output_dir` and `input_dir`. Precedence: file < `-s` < `-o` / `-i`.
- The run writes the configuration it used, defaults included, to `resolved_config.ini`.

`fdtd check run.ini` applies all of these rules without starting MPI.

## Sections

| Section | Sets |
|---|---|
| `[grid]` | domain size, cell size, PML, time steps, structure origin, periodic faces |
| `[parallel]` | MPI ranks per axis |
| `[boxes]` | the total-field / scattered-field box and the two flux boxes |
| `[source]`, `[plane_wave]`, `[wave_superposition]`, `[gaussian_beam]`, `[mode]`, `[dipole.<name>]`, `[wave.<name>]` | the source ([Source Types](../sources/SourceTypes.md)) |
| `[signal.<name>]` | the time signal ([Signal and Time Steps](SignalSettings.md)) |
| `[spectrum]` | the wavelengths every monitor analyses |
| `[geometry]`, `[<shape>]`, `[block.<name>]` | the layers and the structure ([Shapes](../geometry/Shapes.md), [Blocks](../geometry/Blocks.md)) |
| `[material.<name>]`, `[materials] catalog` | the materials by name and where their files are ([Materials](../models/Materials.md)) |
| `[nonlinear]`, `[chi2.<name>]`, `[chi3.<name>]` | second- and third-harmonic generation ([Nonlinear Tensors](../models/NonlinearTensors.md)) |
| `[electrons]`, `[hydrodynamic]` | the free-electron model ([Nonlocal Plasmonics](../models/NonlocalPlasmonics.md), [Hydrodynamics](../models/Hydrodynamics.md)) |
| `[temperature]`, `[two_temperature]`, `[pump_probe]` | the electron temperature ([ITO](../models/ITOmodel.md)) |
| `[output]`, `[multipoles]`, `[probe.<name>]`, `[monitor.<name>]` | what is written ([Viewing the Outputs](Outputs.md)) |
| `[run]`, `[paths]` | wall-clock limit, logging intervals, directories |

## The domain: `[grid]`, `[parallel]`

```ini
[grid]
nx = 120            ; interior cells; pml_x cells are added on each side
ny = 120
nz = 120
dx = 1e-9           ; [m]
dy = 1e-9
dz = 1e-9
pml_x = 20          ; CPML cells per face (default 20)
pml_y = 20
pml_z = 20
time_steps = 8000
```

- `origin_x`, `origin_y`, `origin_z` place the structure, in interior cells; the default `-1` is the
  centre of the domain. Shapes and blocks are positioned from this point, and the layers of
  `[geometry]` meet at `origin_y`.
- `periodic_x`, `periodic_z` make those faces periodic; set the PML of a periodic axis to 0.
- `courant` (default 0.5) fixes the time step, `dt = courant * d / c` with `d` the smallest cell size.
- `staircase = true` assigns one material per cell instead of one per field component.

`[parallel] ranks_x`, `ranks_y`, `ranks_z` split the domain; their product is the number of MPI
ranks the job must run on. `chunks_x = 0 40 90` (start cells) replaces the equal split along an axis.

## The boxes: `[boxes]`

The plane wave is injected on the faces of the total-field / scattered-field (TF/SF) box. Inside it
is the total field, outside only the scattered field. Margins are in cells, counted inward from the
PML:

| Key | Default | Box |
|---|---|---|
| `tfsf_x`, `tfsf_y`, `tfsf_z` | 11 | the TF/SF box |
| `scattering_x`, … | `tfsf - 5` | the scattering-flux box, in the scattered-field region (smaller margin than `tfsf`) |
| `absorption_x`, … | `tfsf + 5` | the absorption-flux box, in the total-field region (larger margin than `tfsf`) |

With periodic x and z only the y margins matter: `tfsf_y` is the injection plane, `scattering_y`
the reflection plane and `absorption_y` the transmission plane.

## The spectrum: `[spectrum]`

```ini
[spectrum]
lambda_min = 300e-9
lambda_max = 800e-9
points = 31
spacing = wavelength        ; or frequency: a uniform step in omega
harmonic_band = none        ; second | third: adds a band at 2w / 3w of every point
```

These keys mean the same thing in every simulation, and every frequency-domain output uses these
points. `points = 1` needs `lambda_min = lambda_max`. To analyse a pump and its second harmonic only:
`lambda_min = lambda_max = <pump>`, `points = 1`, `harmonic_band = second`.

## Layers and structure: `[geometry]`

```ini
[geometry]
first_medium = vacuum       ; where the wave enters, y below the layers
middle_medium = film        ; a layer of middle_thickness that ends at origin_y
middle_thickness = 100e-9
last_medium = glass         ; from origin_y on
shape = sphere              ; the structure; its dimensions are in [sphere]
material = gold_d2cp        ; material of the structure
shift_y = 0                 ; moves the structure along y [m]; positive is toward the source

[sphere]
radius = 20e-9
```

- A material is an integer code: 0 is vacuum, 1 to 99 are dielectrics with index `n<i>` of
  `[materials]`, 100 and above are dispersive materials ([Materials](../models/Materials.md)).
- Every shape has a section named after it. [Shapes](../geometry/Shapes.md) shows each one with its keys.
- `[block.<name>]` sections add boxes of one material (substrates, gratings); they are drawn before
  the shape ([Blocks](../geometry/Blocks.md)).

## Outputs: `[output]`, probes and monitors

```ini
[output]
dft_xy = true               ; frequency-domain E on the plane through the structure
dft_magnetic = true         ; also H on the enabled planes
volume_dft = false          ; E and H in the whole domain (large)
frames_xz = 20              ; a time frame every 20 steps; 0 = off
frames_3d = 0

[probe.gap]                 ; E(t) at one cell (at most 10 probes)
x = 60
y = 60
z = 60

[monitor.near]              ; a box of cells, in the frequency or the time domain
type = frequency            ; frequency | time
center_x = 60
center_y = 60
center_z = 60
size_x = 20
size_y = 20
size_z = 20
step = 1                    ; time monitors: steps between frames
```

The cross-section and reflection / transmission tables are always written.
`[multipoles] exact = true` adds the exact multipole decomposition of the structure, `lwa = true`
the long-wavelength one.

## Run control: `[run]`, `[paths]`

- `[run] wallclock_limit_hours` stops the run cleanly before a time limit
  ([Stopping and Restarting](Restart.md)); 0 means never.
- `[run] trace_interval` is the number of steps between the field trace lines in the log; 0 turns them off.
- `[paths] output_dir` (`-o`) is the output root; `[paths] input_dir` (`-i`) is where data files
  are read from: the MPB files of a mode source, `electron_temperature.h5`.
