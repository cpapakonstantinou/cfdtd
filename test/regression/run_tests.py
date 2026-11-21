#!/usr/bin/env python3
"""Regression tests: run each case with the fdtd driver and compare with EXPECTED_OUTPUT.

Every case directory holds
    run.ini            the configuration
    INIDEF/            input data files, if any (-i)
    EXPECTED_OUTPUT/   reference files, compared value by value:
                         *.h5   HDF5 references
                         *.txt  numeric tables (ART, Q coefficients, ...)

Examples
    python run_tests.py                          # all cases on this machine
    python run_tests.py --testcase ThinFilm --program ../../fdtd
    python run_tests.py --cluster TRILLIUM       # module setup + mpirun on a cluster
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np

HERE = os.path.dirname(os.path.realpath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))

# name: MPI ranks, extra --set overrides, known issue or None
CASES = {
    "EmptyDomain": (8, [], None),
    "EmptyDomain_Gaussian": (8, [], None),
    "ThinFilm": (8, [], None),
    "MieScattering": (8, ["grid.time_steps=1000"], None),
    "EmptyDomainObliqueGaussian": (8, [], None),
    "SubstratePlaneWaveSuperposition": (8, [], None),
    "Anapole": (8, [], None),
}

def load_h5(path):
    """x, y, z stacked; complex components are either one compound dataset (older
    outputs) or the pair <c>_re / <c>_im (current layout, readable by ParaView)."""
    def component(h, c):
        if c in h:
            return h[c][()]
        return h[c + "_re"][()] + 1j * h[c + "_im"][()]
    with h5py.File(path, "r") as h:
        return np.stack([component(h, c) for c in "xyz"], axis=-1)


def load_table(path):
    return np.loadtxt(path, ndmin=2)


def close(ref, new, rtol):
    """Equal within rtol relative to the largest reference magnitude."""
    if ref.shape != new.shape:
        return False
    scale = max(np.nanmax(np.abs(ref)), np.finfo(float).tiny)
    return bool(np.all(np.abs(ref - new) <= rtol * scale))


def compare(case_dir, out_dir, rtol):
    expected = os.path.join(case_dir, "EXPECTED_OUTPUT")
    files = sorted(f for f in glob.glob(expected + "/**/*", recursive=True) if os.path.isfile(f))
    failed = []
    for ref in files:
        rel = os.path.relpath(ref, expected)
        try:
            if rel.endswith(".h5"):
                ok = close(load_h5(ref), load_h5(os.path.join(out_dir, rel)), rtol)
            else:
                ok = close(load_table(ref), load_table(os.path.join(out_dir, rel)), rtol)
        except (OSError, KeyError, ValueError) as e:
            ok = False
            rel += f" ({e.__class__.__name__}: {e})"
        if not ok:
            failed.append(rel)
    print(f"  {len(files) - len(failed)}/{len(files)} reference files match")
    for f in failed[:10]:
        print(f"    differs: {f}")
    return not failed and files


def run_case(name, args, tmp_dir):
    ranks, overrides, issue = CASES[name]
    print(f"+++ {name} ({ranks} ranks) +++")
    if issue and not args.include_known_issues:
        print(f"  skipped, known issue: {issue} (--include-known-issues runs it)")
        return None
    case_dir = os.path.join(HERE, name)
    out_dir = os.path.join(tmp_dir, name)
    cmd = (f"{args.env_setup}{args.mpiexec.format(n=ranks)} {os.path.abspath(args.program)} simulate run.ini "
           f"-o {out_dir} " + ("-i INIDEF " if os.path.isdir(os.path.join(case_dir, "INIDEF")) else "") +
           " ".join(f"--set {o}" for o in overrides))
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "stdout.txt"), "w") as fout:
        rc = subprocess.call(cmd, cwd=case_dir, shell=True, stdout=fout, stderr=subprocess.STDOUT)
    if rc != 0:
        print(f"  simulation failed (exit {rc}), see {out_dir}/stdout.txt")
        return False
    if args.update_references:
        expected = os.path.join(case_dir, "EXPECTED_OUTPUT")
        for ref in glob.glob(expected + "/**/*", recursive=True):
            if os.path.isfile(ref):
                rel = os.path.relpath(ref, expected)
                shutil.copy(os.path.join(out_dir, rel), os.path.join(expected, rel))
        print("  references updated")
        return True
    return bool(compare(case_dir, out_dir, args.tolerance))


MEMORY_CASES = ("ThinFilm", "MieScattering", "Anapole")  # `fdtd memory` against the report of `fdtd mesh`


def array_table(text):
    """{group: (largest rank, total)} of the `memory:` table, and the exact `bytes:` row as "bytes"
    (geometry report or `fdtd memory`)."""
    table = {}
    for key in ("memory", "bytes"):
        rows = [line.split(": ", 1)[1].split(",") for line in text.splitlines() if line.startswith(key + ": ")]
        for row in rows[1:]:  # the first row is the header
            if key == "memory":
                table[row[0]] = (float(row[1]), float(row[2]))
            else:
                table["bytes"] = (int(row[0]), int(row[1]))
    return table


def run_memory_case(name, args, tmp_dir):
    """`fdtd memory` (one process, arrays not allocated) must report exactly the array bytes that the
    real set-up on all ranks (`fdtd mesh`) allocates."""
    ranks, overrides, _ = CASES[name]
    case_dir = os.path.join(HERE, name)
    out_dir = os.path.join(tmp_dir, name + "_mesh")
    sets = " ".join(f"-s {o}" for o in overrides)
    print(f"\n+++ memory estimate: {name} +++")
    mesh = (f"cd {case_dir} && {args.env_setup}{args.mpiexec.format(n=ranks)} {os.path.abspath(args.program)} mesh run.ini -o {out_dir} {sets}"
            f" > {out_dir}.txt 2>&1")
    if subprocess.call(mesh, shell=True, executable="/bin/bash") != 0:
        print(f"  mesh failed, see {out_dir}.txt")
        return False
    est = subprocess.run(f"cd {case_dir} && {args.env_setup}{os.path.abspath(args.program)} memory run.ini {sets}", shell=True, executable="/bin/bash",
                         capture_output=True, text=True)
    if est.returncode != 0:
        print(f"  fdtd memory failed:\n{est.stdout[-600:]}{est.stderr[-600:]}")
        return False
    measured = array_table(open(os.path.join(out_dir, "materials", "geometry_report.txt")).read())
    estimated = array_table(est.stdout)
    good = "bytes" in measured and measured == estimated  # the byte row makes the comparison exact
    if good:
        print(f"  matches the mesh: {measured['all'][1]:.3f} GiB of arrays on {ranks} ranks")
    else:
        for group in sorted(set(measured) | set(estimated)):
            if measured.get(group) != estimated.get(group):
                print(f"  {group}: mesh {measured.get(group)}, memory {estimated.get(group)}")
    return good


CLUSTERS = {
    # name: environment setup, MPI launcher with {n} ranks, scratch variable
    "LOCAL": ("", "mpirun -np {n}", None),
    "LUIS": ("module load iimpi/2020a GCCcore/.10.2.0 hdf5", "srun --mpi=pmi2 --ntasks={n}", "BIGWORK"),
    "HLRN": ("module load gcc/9.3.0 intel/19.1.3 impi/2019.9 hdf5", "srun --ntasks={n}", "WORK"),
    "NIAGARA": ("module load CCEnv; sleep 2; module load arch/avx512; sleep 2; module load StdEnv/2020; "
                "sleep 2; module load intel/2021.2.0 hdf5-mpi", "mpirun -np {n}", "SCRATCH"),
    "NARVAL": ("module load StdEnv/2020 hdf5-mpi", "mpirun -np {n}", "SCRATCH"),
    "TRILLIUM": ("module load StdEnv/2023; sleep 2; module load openmpi hdf5-mpi; sleep 2; module load python",
                 "mpirun -np {n}", "SCRATCH"),
}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--testcase", default="all", help="case name (directory) or 'all'")
    p.add_argument("--program", default=os.path.join(REPO, "fdtd"), help="path to the fdtd executable")
    p.add_argument("--cluster", type=str.upper, default="LOCAL", choices=sorted(CLUSTERS))
    p.add_argument("--mpiexec", help="MPI launcher with {n} for the rank count (overrides --cluster)")
    p.add_argument("--tolerance", type=float, default=1e-10, help="relative to the largest reference value")
    p.add_argument("--include-known-issues", action="store_true", help="also run cases with known issues")
    p.add_argument("--update-references", action="store_true",
                   help="replace each case's EXPECTED_OUTPUT files by this run's outputs")
    p.add_argument("--keeptmpdir", action="store_true", help="keep the simulation outputs")
    args = p.parse_args()

    env, launcher, scratch = CLUSTERS[args.cluster]
    args.env_setup = env + "; " if env else ""
    args.mpiexec = args.mpiexec or launcher
    base = os.environ.get(scratch) if scratch else None
    tmp_dir = tempfile.mkdtemp(prefix="fdtd_regression_", dir=base)
    print(f"running in {tmp_dir}")

    names = list(CASES) if args.testcase == "all" else [c for c in CASES if c == args.testcase]
    if not names:
        p.error(f"unknown test case '{args.testcase}' (have: {', '.join(CASES)})")
    results = {n: run_case(n, args, tmp_dir) for n in names}
    if not args.update_references:
        for n in names:
            if n in MEMORY_CASES:
                results["memory: " + n] = run_memory_case(n, args, tmp_dir)

    if not args.keeptmpdir:
        shutil.rmtree(tmp_dir)
    ran = {n: r for n, r in results.items() if r is not None}
    print(f"\n{sum(ran.values())}/{len(ran)} cases passed" + (f", {len(results) - len(ran)} skipped" if len(ran) < len(results) else ""))
    sys.exit(0 if all(ran.values()) else 1)


if __name__ == "__main__":
    main()
