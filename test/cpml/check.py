#!/usr/bin/env python3
"""CPML reflection: max |E_test - E_reference| / peak at probes 4 cells from the CPML.

  cd test/cpml
  mpirun -n 2 ../../fdtd simulate domain_24.ini  -o out_24
  mpirun -n 2 ../../fdtd simulate domain_160.ini -o out_160
  python3 check.py out_24 out_160

2026-11-06: face -79.4 dB, corner -66.6 dB (-66.7 after the alpha grading fix); identical on 1 rank, 2 ranks and with a
rank boundary inside the CPML (parallel.chunks_x = "0 5").
"""
import sys, h5py, numpy as np
test, ref = sys.argv[1], sys.argv[2]
def series(d, name):
    out = []
    t = 1
    while True:
        try:
            out.append(h5py.File(f'{d}/frames/E_Probe_{name}_Frame_{t}.h5')['z'][()].ravel()[0])
        except OSError:
            return np.array(out)
        t += 1
for name in ['face', 'corner']:
    a, b = series(test, name), series(ref, name)
    e = np.abs(a - b).max() / np.abs(b).max()
    print(f'{name:7s} {20 * np.log10(e):6.1f} dB')
