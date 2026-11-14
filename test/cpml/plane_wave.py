#!/usr/bin/env python3
"""CPML reflection at normal incidence: a sheet of z-dipoles in a domain periodic in x and z
launches plane waves along +-y; a probe 4 cells from the far y CPML is compared with the same
probe in a domain long enough that no CPML reflection reaches it in the time window.

  python3 plane_wave.py write [pml] [n]   # writes test.ini and reference.ini, background index n
  mpirun -n 1 ../../fdtd simulate test.ini -o out_test
  mpirun -n 1 ../../fdtd simulate reference.ini -o out_ref
  python3 plane_wave.py check out_test out_ref

5 fs pulse at 800 nm, 20 nm cells, 1600 steps. 2026-11-06 (after the alpha grading fix):
  vacuum:  4 cells -48.3 dB, 8 cells -82.8 dB, 16 cells -100.9 dB
  n = 1.5: 8 cells -79.1 dB (-83.7 dB with sigma / n; vacuum at the same cells per wavelength: -84.2)
  n = 3.5: 8 cells -60.3 dB (-59.2 with sigma / n), 16 cells -90.1 dB; vacuum at the same cells per
           wavelength (70 nm cells): -54.0 dB. The loss in dense media is resolution, not the grading.
The alpha grading fix (integer m / N at the E positions) changed none of these by more than 0.1 dB.
"""
import sys

def ini(ny, y_sheet, y_probe, pml, n):
    medium = '' if n == 1 else f'first_medium = 4\nmiddle_medium = 4\nlast_medium = 4\n[materials]\nn4 = {n}\n'
    s = f"""[grid]
nx = 4
ny = {ny}
nz = 4
dx = 20e-9
dy = 20e-9
dz = 20e-9
pml_x = 8
pml_y = {pml}
pml_z = 8
periodic_x = true
periodic_z = true
time_steps = 1600
[boxes]
tfsf_x = 0
tfsf_y = 2
tfsf_z = 0
scattering_x = 0
scattering_y = 1
scattering_z = 0
absorption_x = 0
absorption_y = 1
absorption_z = 0
[source]
kind = dipole
pulse = gaussian
center_wavelength = 800e-9
pulse_width = 5e-15
[spectrum]
lambda_min = 800e-9
lambda_max = 800e-9
points = 1
[geometry]
shape = none
{medium}"""
    for name, y in (('far', y_probe), ('inc', y_sheet + 4)):
        s += f'[monitor.{name}]\ntype = time\ncenter_x = 2\ncenter_y = {y}\ncenter_z = 2\nsize_x = 1\nsize_y = 1\nsize_z = 1\n'
    for x in range(4):
        for z in range(4):
            s += f'[dipole.d{x}{z}]\nx = {x}\ny = {y_sheet}\nz = {z}\ne_x = 0\ne_y = 0\ne_z = 1\n'
    return s

def series(d, name):
    import h5py, numpy as np
    out, t = [], 1
    while True:
        try:
            out.append(h5py.File(f'{d}/frames/E_Probe_{name}_Frame_{t}.h5')['z'][()].ravel()[0])
        except OSError:
            return np.array(out)
        t += 1

if sys.argv[1] == 'write':
    pml = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    n = float(sys.argv[3]) if len(sys.argv) > 3 else 1
    open('test.ini', 'w').write(ini(40, 20, 36, pml, n))
    open('reference.ini', 'w').write(ini(1800, 900, 916, pml, n))
else:
    import numpy as np
    a, b = series(sys.argv[2], 'far'), series(sys.argv[3], 'far')
    print(f'reflection {20 * np.log10(np.abs(a - b).max() / np.abs(b).max()):.1f} dB')
