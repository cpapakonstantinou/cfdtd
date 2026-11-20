# 3D-FDTD

A parallel (MPI) finite-difference time-domain solver for nano-optics: scattering and absorption
of nanoantennas, layered and periodic structures, dispersive and nonlinear materials.

- **Materials:** dielectrics, Drude, Drude + critical points, Lorentz, modified Lorentz, magnetic
  Lorentz; nonlocal and hydrodynamic free electrons; a two-temperature model; χ² and χ³ tensors.
- **Sources:** plane wave (TF/SF), superposed plane waves through layers, focused Gaussian beam,
  point dipoles, waveguide modes from MPB.
- **Monitors:** scattering / absorption cross-sections, reflection / transmission, DFT planes,
  boxes and volumes, time frames and probes, near-to-far faces, multipoles.

The user guide is in `docs/` (start at `docs/index.md`); `mkdocs build` renders it as a site.

## Build

Linux, an MPI C++23 compiler (`mpicxx`) and parallel HDF5. On Debian / Ubuntu:
`libhdf5-openmpi-dev`. HDF5 is found through `pkg-config hdf5`; if that name is the serial library
on your machine the build stops with an error, and `make HDF5_PKG=hdf5-openmpi` names the parallel one.

```sh
make                      # ./fdtd, optimised for this machine (-march=native)
make MARCH= OPT=-O2       # portable build; results are bit-comparable across machines
```

## Run

A run is described by one INI file.

```sh
./fdtd template > run.ini                    # every key, with its default and a description
./fdtd check run.ini                         # validate; no MPI
./fdtd memory run.ini                        # memory per rank and in total, without allocating it
mpirun -n 8 ./fdtd mesh     run.ini -o out   # build the mesh only: material maps and a report
mpirun -n 8 ./fdtd simulate run.ini -o out   # run
```

- `-o DIR` sets the output directory, `-i DIR` the directory of data files (mode profiles, an
  electron temperature), and `-s section.key=value` overrides any key.
  Precedence: file < `-s` < `-o` / `-i`.
- The number of MPI ranks must equal `ranks_x * ranks_y * ranks_z` of `[parallel]`.
- An unknown or misspelled key is an error. The run writes the configuration it used to
  `resolved_config.ini`.
- Structured outputs are HDF5 files with an XDMF companion for ParaView; tables are text.

`examples/` holds complete configurations; `test/regression/` holds small ones with reference
results. See `docs/tutorial/CodeTutorial.md` for a first run.

## Test

```sh
make tests check-tests MPIEXEC="mpirun -n 3"    # unit tests
make regression                                  # 7 small runs against stored references
make check-layers                                # module dependency rules
```

`docs/dev/Tests.md` describes them.

## Documentation

```sh
pip install mkdocs python-markdown-math
mkdocs build              # user guide -> mkdocs_output/
make doxygen              # API reference -> docs/output/ (generates the Doxyfile)
```

## Layout

| Path | Content |
|---|---|
| `inc/<module>/`, `src/<module>/` | the solver, one directory per module (`docs/dev/Architecture.md`) |
| `tools/` | the layer check, the migration script for old input files |
| `test/` | unit tests, regression cases, CPML and restart checks |
| `examples/` | complete configurations |
| `docs/` | user guide and developer notes; `docs/dev/Roadmap.md` records the state of the code and every intentional change of results |
