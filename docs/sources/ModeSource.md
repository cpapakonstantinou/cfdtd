# Mode Source

`[source] type = mode` injects a guided mode on a plane normal to x, y or z, travelling one way. The
mode comes from [MPB](https://mpb.readthedocs.io): its transverse E and H on the waveguide's
cross-section, written by MPB as HDF5 files. No other tool is involved.

## 1. Solve the mode in MPB

Solve the cross-section (a 2D cell, no size along the waveguide) at the run's wavelength with
`find_k`, and write E and H of the band:

```python
ms = mpb.ModeSolver(geometry_lattice=mp.Lattice(size=mp.Vector3(cell_x, cell_y)),
                    geometry=[...], default_material=..., resolution=50, num_bands=1,
                    filename_prefix="strip")
omega = 1 / wavelength                       # in units of a
k = ms.find_k(mp.NO_PARITY, omega, 1, 1, mp.Vector3(0, 0, 1), 1e-7,
              omega * n_core, omega * n_clad, omega * n_core,
              mpb.fix_efield_phase, mpb.output_efield, mpb.fix_hfield_phase, mpb.output_hfield)
print("effective_index =", k[0] / omega)
```

This writes `strip-e.k01.b01.h5` and `strip-h.k01.b01.h5`. `examples/mode_source/strip_mode.py`
and `examples/grating_antenna/waveguide_mode.py` are complete inputs.

- Make the cell large enough that the field is negligible at its edges: the profile is set to zero
  outside it.
- MPB's resolution is independent of the FDTD grid: the profile is interpolated (bilinearly) onto the
  Yee positions of the injection plane. A finer MPB grid than the FDTD grid is fine.
- Use the materials the FDTD run uses at that wavelength (the same refractive indices).

## 2. Configure the run

```ini
[source]
type = mode
signal = gaussian              ; a signal with a carrier at the mode's wavelength

[signal.gaussian]
center_wavelength = 1.55e-6
width = 8e-15

[mode]
efield = strip-e.k01.b01.h5    ; relative to [paths] input_dir (-i)
hfield = strip-h.k01.b01.h5
unit_length = 1e-6             ; MPB's a in metres
effective_index = 1.652        ; k / omega from find_k
propagation_axis = x           ; the mode travels along +x
position = -2.4e-6             ; injection plane, from the structure centre
; mpb_x_axis = z               ; FDTD axes of MPB's x and y; defaults below
; mpb_y_axis = y
; offset_x = 0                 ; MPB cell centre from the structure centre [m]
```

| `propagation_axis` | MPB x along (default) | MPB y along (default) |
|---|---|---|
| `x` | z | y |
| `y` | x | z |
| `z` | x | y |

The defaults assume MPB's y is "up" for a waveguide along x or z, with the layers stacked along y
as the antenna shapes expect. MPB's cell is centred on the structure centre (`[grid] origin_*`) plus
`offset_*`; MPB's grid points are taken at `-L/2 + i L/N` from its cell centre. The run logs where
the mode's |E|² centroid lands, as a check of the placement:

```
sources: mode source: |E_t|^2 centroid at MPB (x, y) = (0.0118, -0.0009) um from the cell centre, i.e. z +0.0118 um, y -0.0009 um from the structure centre
```

The structure the mode travels in is drawn with `[block.<name>]` sections (boxes of one material,
in file order, before the `[geometry]` shape). A block whose bounds are left out reaches past the
domain on that side; with the mode source, the material maps continue through the PML, so a
waveguide that leaves the domain is absorbed rather than reflected.

## 3. How it is injected

The plane is a total-field / scattered-field boundary: downstream of it the total field, upstream
only the scattered field. With (n, u, v) cyclic (n the propagation axis) and p the E plane:

```
E_u(p)       += dt / (eps dn)  H_v,inc(p - 1/2)       H_v(p - 1/2) += dt / (mu0 dn) E_u,inc(p)
E_v(p)       -= dt / (eps dn)  H_u,inc(p - 1/2)       H_u(p - 1/2) -= dt / (mu0 dn) E_v,inc(p)
```

The incident fields are the mode profile times the signal; the H plane leads the E plane by half a
cell at `effective_index`, and E is taken half a time step before the H the E update uses. eps is
the material of each cell of the plane, so the correction is right in the core and in the cladding.

- **Amplitude:** the profile's peak |E_t| is `[source] peak_field` times the signal envelope. H is
  MPB's H divided by the vacuum impedance (MPB's units).
- **Phase and direction:** MPB's fields are complex with an arbitrary phase. E and H are each turned
  real by removing the phase of their largest sample. The sign of H is then chosen so that the power
  flows along +n. This makes the input independent of the direction MPB solved for and of how its
  axes are mapped.
- **One way:** with the profile computed on a 20 nm MPB grid and injected on a 40 nm FDTD grid (Si3N4
  strip in SiO2), the field upstream of the plane is -46 dB of the forward field. The remainder is the
  difference between MPB's mode and the FDTD grid's own mode. A finer FDTD grid lowers it.

## 4. Limits

- One profile at one wavelength. A short pulse injects that profile at every frequency of its band;
  the mode changes slowly with wavelength in most waveguides, but keep the band narrow when the
  profile matters.
- Lossless dielectric materials on the injection plane. A Drude-type metal on the plane uses its
  instantaneous coefficient; other dispersive models get no E correction.
- No restart support yet (`-r`).

## Worked examples

**`examples/mode_source/waveguide.ini`**: the TE mode of a 1 µm × 400 nm Si3N4 strip in SiO2,
injected along +x. Run `strip_mode.py` with MPB, copy the printed effective index into the file, then

```
mpirun -n 4 ./fdtd simulate examples/mode_source/waveguide.ini -i examples/mode_source -o out/waveguide
```

`fields/E_1550nm_XY.h5` (side view through the strip) shows the mode downstream of x = -2.4 µm and
nothing upstream.

**`examples/grating_antenna/grating_antenna.ini`**: a gold dipole fed by a Si3N4 grating coupler. A
3 µm box on a Si handle, a 400 nm Si3N4 input waveguide that ends in a fully etched grating (8 teeth,
1.15 µm period), SiO2 cladding of adjustable height, and the dipole on the cladding above the third
period. The waveguide's TE mode is injected along +x; the grating sends it up at about 10° in SiO2
onto the antenna. The header of the file lists the two keys that set the cladding height, and the
cost of the grid (42 M cells at 20 nm, 340 M at 10 nm).
