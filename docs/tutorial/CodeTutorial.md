# Running a Simulation

This page builds the solver and makes two small runs: an empty domain and a thin film. Both
configurations are in `test/regression/` and finish in under a minute on a laptop.

## Build

The solver needs Linux, an MPI C++23 compiler (`mpicxx`) and parallel HDF5 (Debian / Ubuntu:
`libhdf5-openmpi-dev`).

```sh
make            # builds ./fdtd
```

`./fdtd --help` lists the commands:

| Command | Does |
|---|---|
| `fdtd template` | prints a configuration file with every key, its default and a description |
| `fdtd check run.ini` | validates the file and tests that the time steps cover the signal; no MPI |
| `fdtd memory run.ini` | counts the memory the run needs, per rank and in total, without allocating it; one process |
| `fdtd mesh run.ini -o out` | builds the mesh only and writes the material maps and a geometry report |
| `fdtd simulate run.ini -o out` | builds the mesh and runs |

## First run: an empty domain

`test/regression/EmptyDomain/run.ini` is a plane wave through vacuum:

```ini
[grid]
nx = 80                 ; interior cells; the PML is added on top
ny = 120
nz = 80
dx = 1e-9               ; cell size [m]
dy = 1e-9
dz = 1e-9
pml_x = 20              ; CPML cells per face
pml_y = 20
pml_z = 20
time_steps = 500

[parallel]
ranks_x = 2             ; 2 x 2 x 2 = 8 MPI ranks
ranks_y = 2
ranks_z = 2

[boxes]
tfsf_x = 20             ; total-field / scattered-field box: cells from the PML
tfsf_y = 20
tfsf_z = 20

[spectrum]
lambda_min = 300e-9     ; wavelengths analysed by every monitor
lambda_max = 800e-9
points = 31

[output]
frames_xy = 10          ; an E-field frame of the XY plane every 10 steps
```

Every key that is not in the file has its default: the source is a plane wave travelling along +y,
polarised along z, with a broadband raised-cosine pulse.

Validate the file, then run it with as many ranks as `[parallel]` asks for:

```sh
./fdtd check test/regression/EmptyDomain/run.ini
mpirun -n 8 ./fdtd simulate test/regression/EmptyDomain/run.ini -o out/empty
```

The number of ranks given to `mpirun` must equal `ranks_x * ranks_y * ranks_z`; the run stops with
a message otherwise. On a machine with fewer than 8 cores add `--oversubscribe` (Open MPI), or
lower the rank counts: `-s parallel.ranks_x=1 -s parallel.ranks_z=1` gives 2 ranks.

The output directory now holds:

| Path | Content |
|---|---|
| `resolved_config.ini` | the configuration the run used, with every default filled in |
| `reports/art_tfsf_500.txt` | transmission, reflection and absorption per wavelength |
| `reports/cross_sections.txt` | scattering, absorption and extinction efficiencies |
| `reports/source_signal.txt` | the source signal |
| `frames/Efield_XY_<n>.h5` (+ `.xmf`) | the time frames |
| `materials/MaterialIndex_*.h5` | the material maps |

In an empty domain the transmission is 1 and the reflection is 0 (about 1e-35 here). Open
`frames/Efield_XY_00010.xmf` in ParaView to see the pulse; ParaView loads the numbered files as a
time series. [Viewing the Outputs](Outputs.md) describes the files.

## Second run: a thin film

`test/regression/ThinFilm/run.ini` is a 100 nm dielectric slab, periodic in x and z:

```ini
[grid]
pml_x = 0               ; no PML on the periodic faces
pml_z = 0
periodic_x = true
periodic_z = true

[geometry]
shape = none            ; no structure, only the layers
middle_medium = film    ; the material film fills the middle layer
middle_thickness = 100e-9

[material.film]
model = dielectric
n = 2.46                ; refractive index of the film
```

The domain along y is three layers: the first medium where the wave enters, the middle layer of
`middle_thickness` ending at the structure origin, and the last medium. A material is named and
defined in a `[material.<name>]` section or a catalog file; `vacuum` is built in and is the default.

```sh
mpirun -n 8 ./fdtd simulate test/regression/ThinFilm/run.ini -o out/film
python3 test/regression/ThinFilm/check.py out/film
```

`check.py` plots the transmission and reflection of `reports/art_4000.txt` against the
analytical thin-film result and prints the deviations (about 0.2 % in T).

To explore: change `middle_thickness` or the film's `n` from the command line, without editing the file,

```sh
mpirun -n 8 ./fdtd simulate test/regression/ThinFilm/run.ini -o out/film2 -s material.film.n=1.5
```

## Next steps

- A nanoparticle: `test/regression/MieScattering/run.ini` (a 20 nm gold sphere, `shape = sphere`,
  `material = gold_d2cp`). [Shapes](../geometry/Shapes.md) lists the structures and their keys.
- Look at the mesh before running: `mpirun -n 8 ./fdtd mesh run.ini -o out` writes
  `materials/MaterialIndex_*.h5` and `materials/geometry_report.txt` in a few seconds.
- Choose the number of time steps: [Signal and Time Steps](SignalSettings.md).
- Larger examples are in `examples/`.

## How much memory will it need?

`fdtd memory` answers that before the job is submitted. It is one process, needs no `mpirun`, and
runs on a laptop for a job of any size:

```sh
./fdtd memory run.ini
```

```
layout: ranks_x,ranks_y,ranks_z,ranks
layout: 4,4,4,64
memory: group,largest_rank_gib,total_gib
memory: fields,1.116,71.440
memory: material_maps,0.187,11.953
memory: material_models,1.639,26.542
memory: cpml,0.054,1.741
memory: source,0.000,0.011
memory: monitors,0.234,9.169
memory: other,0.000,0.000
memory: all,3.022,120.857
bytes: largest_rank,total
bytes: 3244590976,129769474048
job: smallest_rank_gib,overhead_per_rank_gib,largest_rank_gib,total_gib
job: 1.628,0.059,3.081,124.633
```

A 1000³ domain with a 400 nm gold sphere; the command took 11 s and 0.28 GiB. 1 GiB = 1024³ bytes,
the unit of a scheduler's `--mem=...G`.

- **Largest rank** is what the scheduler's memory per task has to cover; **total** is the job.
  Ranks are uneven: the ones that hold the structure, PML or monitor planes need more.
- **Trying a layout** costs nothing: `./fdtd memory run.ini -s parallel.ranks_x=8 -s parallel.ranks_y=4`.
  More ranks add the process overhead and one layer of neighbour cells per internal face, so the
  total grows slowly with the rank count while the largest rank shrinks.
- **The numbers are counted, not modelled.** The command runs the solver's own set-up for every rank
  in turn, up to and including the mesh, with arrays that take no memory until they are written.
  The array figures are therefore the ones a real run allocates; `fdtd mesh` prints the same table
  in `materials/geometry_report.txt`, and a run logs its peak in the `arrays:` table when it ends.
- **Not counted:** the process overhead, MPI, HDF5 and the other libraries, is a measured
  constant, 0.04 to 0.065 GiB per rank with Open MPI; output buffers at the end of a run come on
  top, about 5 % in a run with a volume DFT.

## On a cluster

Run in the scratch file system and start the solver with the scheduler's launcher. With SLURM:

```sh
#!/bin/bash
#SBATCH --ntasks=64
#SBATCH --time=12:00:00
#SBATCH --signal=B:USR1@300        # lets the run save its state before the limit

srun ./fdtd simulate run.ini -o $SCRATCH/run42
```

`--ntasks` must equal the product of the three rank counts. [Stopping and Restarting](Restart.md)
explains the signal line.
