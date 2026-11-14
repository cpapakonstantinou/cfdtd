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
| `reports/` | text tables: `Qcoeff_*`, `ART-TFSF_*`, `ART_*`, `NL_Efficiency_*`, `Signal.txt`, `DFT_module_phase.txt`, the near/far-field point files, the 1D profiles (`1D_Ex.txt` ...), and the logs |
| `fields/` | frequency-domain data: DFT planes, boxes and volumes, the near-to-far faces |
| `frames/` | time-domain data: movie frames `Efield_<plane>_<number>.h5` (the number last, so ParaView loads them as a series) and the probe traces |
| `materials/` | the material index maps (`MaterialIndex_*`), the geometry report |
| `state/` | application state dumps (temperatures; the restart file later) |

A directory the run does not write to is not kept. `resolved_config.ini` is in the root.

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
