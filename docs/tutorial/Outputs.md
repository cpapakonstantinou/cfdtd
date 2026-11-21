# Viewing the Outputs

Every structured output (DFT planes and volumes, time frames, material masks, custom
monitor boxes, near-to-far faces, temperature) is an HDF5 file. Next to each
`<name>.h5` the solver writes `<name>.xmf`, an XDMF description of the same data, so
the file opens directly in ParaView or VisIt. Turn this off with

```ini
[output]
xdmf = false
```

## The output tree

| Directory | Content |
|---|---|
| `reports/` | text tables, see [Reports](#reports), and `run.log` |
| `fields/` | frequency-domain data: DFT planes, boxes and volumes, the near-to-far faces |
| `frames/` | time-domain data: movie frames `Efield_<plane>_<number>.h5` (the number last, so ParaView loads them as a series) and the probe traces |
| `materials/` | the material index maps (`MaterialIndex_*`), the geometry report |
| `state/` | application state dumps (temperatures; the restart file later) |

A directory the run does not write to is not kept. `resolved_config.ini` is in the root.

## Reports

Every text report starts with comment lines: what the run was, then the column header.

```
# 3D-FDTD 9e76cd2, double precision
# config: run.ini
# override: grid.time_steps=1000
# source: plane_wave, gaussian, polarization z, 0 deg toward the other transverse axis, propagation y
# spectrum: 690 to 785 nm, 5 points, frequency spacing, harmonic band third
# grid: 500 x 500 x 500 cells of 1 x 1 x 1 nm, cpml 16 16 16
# time: dt 1.667820e-18 s, 65000 steps
# LAMBDA[nm]	 Qscatt	 Qabs	 Qextin
7.850000e+02	 ...
```

`numpy.loadtxt(path)` and gnuplot skip them; with pandas,
`pandas.read_csv(path, comment="#", sep=r"\s+", header=None)`.

| File | Content |
|---|---|
| `cross_sections.txt` | scattering, absorption and extinction efficiencies per wavelength |
| `cross_sections_harmonic.txt` | the same at the harmonic, hydrodynamic runs |
| `cross_sections_sigma.txt` | with the absorption from sigma \|E\|^2, `[output] sigma_absorption` |
| `art_tfsf_<n>.txt` | T, R, A from the flux boxes at step n; `art_tfsf_harmonic_<n>.txt`, `art_tfsf_sigma.txt` |
| `art_<n>.txt`, `art_sigma.txt` | T, R, A from the planes of periodic runs |
| `harmonic_efficiency.txt` | harmonic conversion efficiency from the flux boxes; `harmonic_efficiency_art.txt` from the ART-TFSF estimate, `harmonic_efficiency_periodic.txt` from the planes |
| `multipoles_<method>.txt` | multipole moments, `exact` or `lwa`; `multipoles_<method>_full_current.txt` at the harmonic |
| `multipoles_<method>_cross_sections.txt` | the cross section of each multipole |
| `source_signal.txt` | the source signal and its carrier, one row per step |
| `source_spectrum.txt` | the source spectrum the DFTs are normalised by |
| `far_field_forward.txt`, `far_field_backward.txt`, `near_field.txt` | Ex and Ez, amplitude and phase, at one point |
| `aux_grid_ex.txt` ... `aux_grid_hz.txt` | the fields of the 1D auxiliary grid |
| `run.log` | the end of the run, with the elapsed time |

The probe traces are `frames/probe_<x>_<y>_<z>.txt`, the tracker spectra
`frames/tracker_<n>.txt`, with the same comment lines.

## Layout of the HDF5 files

| Dataset | Content |
|---|---|
| `x`, `y`, `z` | real field components (time frames, masks) |
| `x_re`, `x_im`, ... | real and imaginary part of a complex component (DFT outputs) |
| `Ex_re`, `Hy_im`, ... | the tangential fields on a near-to-far face |

File attributes: `domain_size` (the physical domain in cells), `plane` (`XY`, `XZ`, `YZ`)
or `face` (`x`, `y`, `z`), `wavelength` [m]. In Python:

```python
import h5py
with h5py.File("fields/E_500nm_XZ.h5") as f:
    ex = f["x_re"][()] + 1j * f["x_im"][()]   # indexed [x][z]
```

## ParaView

Open the `.xmf` file (reader: *XDMF Reader*) and press *Apply*; every dataset is a point
array of a uniform grid with the cell size of the run. Combine the real and imaginary
parts with the *Calculator* filter, for example `sqrt(x_re^2 + x_im^2)` for |E_x|.

The viewer reverses the index order of the HDF5 arrays (the last index runs along its X
axis). For volumes, ParaView's X, Y, Z are the solver's z, y, x; for planes, ParaView's
Y and Z are the plane's second and first axis (for an `XZ` plane: Y = z, Z = x). Each
`.xmf` states its mapping in a comment.
