#!/usr/bin/env python3
"""Compare two output trees of `fdtd simulate` bit for bit.

    test/compare_runs.py REFERENCE NEW

Every file of REFERENCE must exist in NEW with the same content:
  *.h5   every dataset equal (NaN equals NaN); a missing .xmf companion is counted
  text   the numeric rows equal; files without numbers are compared byte by byte and
         reported separately ("text-only differ")
Logs (*.log, *.out) and resolved_config.ini are skipped: they carry the version and timings.
Exit status 0 when nothing differs.
"""
import glob
import os
import sys

import h5py
import numpy as np


def numbers(path):
    rows = []
    for line in open(path, errors="replace"):
        try:
            rows.append([float(x) for x in line.split()])
        except ValueError:
            pass  # a header or a text line
    return rows


def same_rows(u, v):
    return len(u) == len(v) and all(np.array_equal(np.array(p), np.array(q), equal_nan=True) for p, q in zip(u, v))


def main(ref, new):
    compared, differ, text, no_xmf = 0, [], [], []
    for f in sorted(glob.glob(ref + "/**/*", recursive=True)):
        rel = os.path.relpath(f, ref)
        if not os.path.isfile(f) or rel.endswith(("resolved_config.ini", ".log", ".out", ".xmf")):
            continue
        g = os.path.join(new, rel)
        if not os.path.exists(g):
            differ.append("missing " + rel)
            continue
        compared += 1
        if f.endswith(".h5"):
            with h5py.File(f) as x, h5py.File(g) as y:
                for k in x:
                    if k not in y:
                        differ.append(f"{rel}:{k} missing dataset")
                        continue
                    u, v = x[k][()], y[k][()]
                    if u.shape != v.shape or not np.array_equal(u, v, equal_nan=True):
                        differ.append(f"{rel}:{k} differs")
            if os.path.exists(f[:-3] + ".xmf") and not os.path.exists(g[:-3] + ".xmf"):
                no_xmf.append(rel)
        else:
            u, v = numbers(f), numbers(g)
            if u or v:
                if not same_rows(u, v):
                    differ.append(f"{rel} numbers differ ({len(u)} vs {len(v)} rows)")
            elif open(f, "rb").read() != open(g, "rb").read():
                text.append(rel)
    print(f"{compared} files compared, {len(differ)} differ, {len(text)} text-only differ, {len(no_xmf)} h5 without xmf")
    for line in differ[:12]:
        print("  ", line)
    for name in sorted(set(os.path.basename(t) for t in text))[:6]:
        print("   text:", name)
    return 1 if differ else 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
