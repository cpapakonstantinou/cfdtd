# Testing

Four kinds of checks. Run the first three before every commit; they take a few minutes.

```sh
make MARCH= OPT=-O2                               # the portable build the references were made with
make tests check-tests MPIEXEC="mpirun -n 3"      # unit tests
make regression                                   # 7 runs against stored references, 3 memory checks
make check-layers                                 # module dependency rules
```

On a machine with fewer cores than a case has ranks, pass the launcher with Open MPI's
`--oversubscribe`: `make regression REGRESSION_ARGS='--mpiexec "mpirun --oversubscribe -np {n}"'`.

## Unit tests

`test/test_*.cc`, built by `make tests` and run by `make check-tests`:

| Test | Covers |
|---|---|
| `test_config` | the INI reader, every validation rule, command-line overrides, the template round trip, shape sections, the χ tensors |
| `test_machine` | the state machine of the driver: the order of the phases, error and restart paths |
| `test_carray` | the array and view types |
| `test_geometry` | the set operations and the rendering of the geometry layers; the anapole against its point-wise definition |
| `test_h5io` | parallel HDF5 output, on 1 rank and under `MPIEXEC` |

## Regression cases

`test/regression/run_tests.py` runs each case with `fdtd simulate` and compares the files in its
`EXPECTED_OUTPUT/` value by value (relative tolerance 1e-10). For three cases it also checks that
`fdtd memory` reports exactly the array bytes that `fdtd mesh` allocates on all ranks.

| Case | What it is |
|---|---|
| `EmptyDomain` | a plane wave through vacuum: T = 1, R = 0 |
| `EmptyDomain_Gaussian` | the same with a Gaussian pulse |
| `EmptyDomainObliqueGaussian` | two oblique plane waves (`wave_superposition`) |
| `ThinFilm` | a 100 nm dielectric slab, periodic in x and z |
| `MieScattering` | a 20 nm gold sphere (Drude + critical points) |
| `SubstratePlaneWaveSuperposition` | superposed waves through a layered medium (the 1D simulations) |
| `Anapole` | the anapole nanoantenna on glass (a shape of the geometry layers, layered TF/SF) |

```sh
cd test/regression
python3 run_tests.py                         # all cases
python3 run_tests.py --testcase ThinFilm     # one case
python3 run_tests.py --keeptmpdir            # keep the run directories for inspection
python3 run_tests.py --update-references     # replace the references with the new results
```

Three cases have a `check.py` that compares a run with the analytical result:
`python3 <case>/check.py <output directory>` for `EmptyDomain`, `ThinFilm` and `MieScattering`
(run the Mie case with the 8000 steps of its `run.ini`; the suite shortens it to 1000). The numbers
to expect are in the [Git Guide](GitGuide.md). `EmptyDomain_Gaussian` and the two superposition
cases are regression-only: their pulses are longer than the runs, so their spectra are not physical.

The references are tied to the build flags: `-march=native` uses fused multiply-add and gives results
that differ in the last digits. Use `make MARCH= OPT=-O2`, and `make clean` when switching.

**Updating references.** A refactoring must not change them. When a change is meant to move
results (a bug fix, a corrected formula), find out how far they move and why, record it in
`Roadmap.md` §7, and only then run `--update-references`.

## Comparing two builds bit for bit

A refactoring is checked by running the same configurations with the build before and after, and
comparing the complete output trees:

```sh
mpirun -n 2 ./fdtd_before simulate case.ini -o out/before
mpirun -n 2 ./fdtd        simulate case.ini -o out/after
python3 test/compare_runs.py out/before out/after
```

`compare_runs.py` compares every dataset of every HDF5 file and every numeric table exactly and
prints the files that differ. Cover what the change touches: the material models it can reach, a
periodic and a non-periodic case, more than one rank layout, and a nonlinear case when the update
path is involved. `git worktree add ../before <commit>` gives the build to compare with.

## Other checks

| Script | Checks |
|---|---|
| `test/restart_test.sh` | an interrupted and restarted run against an uninterrupted one, bit for bit |
| `test/edge_layouts.sh` | rank boundaries exactly on the CPML interfaces against the 1-rank run |
| `test/cpml/check.py`, `test/cpml/plane_wave.py` | CPML reflection in dB (dipole near a corner; plane wave at normal incidence); each file states the expected values |
| `tools/check_layers.py --files` | lists every dependency that goes against the module order |
| `./fdtd mesh` + `docs/geometry/make_shape_figures.py` | the mesh of every shape, with its keys drawn on it from their values: a dimension line off an edge means a key does not do what the figure says |

## A warning-free build

The build uses `-Wall -Wextra`. `make clean && make MARCH= OPT=-O2` must print no warning.
