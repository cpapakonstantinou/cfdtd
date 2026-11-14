"""MPB input: the fundamental TE mode of the Si3N4 input waveguide in SiO2 at 1550 nm, for examples/grating_antenna (3 um wide, 400 nm thick).

    python waveguide_mode.py        (needs MPB's Python interface: conda install -c conda-forge pymeep)

Writes waveguide-e.k01.b01.h5 and waveguide-h.k01.b01.h5 (the mode's E and H on the cross-section) and
prints the effective index for [mode] effective_index. Lengths are in units of a = 1 um, so
[mode] unit_length = 1e-6. MPB's x is the strip's width and y its height; with the injection
along x in the FDTD run they lie along z and y ([mode] mpb_x_axis = z, mpb_y_axis = y, the defaults).
"""
import meep as mp
from meep import mpb

wavelength = 1.55          # um
width, height = 3.0, 0.4   # input waveguide, um
cell_x, cell_y = 5.0, 2.4  # cross-section, um: large enough for the field to vanish at the edges
n_core, n_clad = 2.0, 1.444

ms = mpb.ModeSolver(
    geometry_lattice=mp.Lattice(size=mp.Vector3(cell_x, cell_y)),   # no size along z: the waveguide axis
    geometry=[mp.Block(size=mp.Vector3(width, height, mp.inf), material=mp.Medium(index=n_core))],
    default_material=mp.Medium(index=n_clad),
    resolution=50,         # 20 nm; the FDTD run interpolates onto its own grid
    num_bands=1,
    filename_prefix="waveguide",
)
omega = 1 / wavelength
k = ms.find_k(mp.NO_PARITY, omega, 1, 1, mp.Vector3(0, 0, 1), 1e-7,
              omega * n_core, omega * n_clad, omega * n_core,
              mpb.fix_efield_phase, mpb.output_efield, mpb.fix_hfield_phase, mpb.output_hfield)
print(f"effective_index = {k[0] / omega:.6f}")
