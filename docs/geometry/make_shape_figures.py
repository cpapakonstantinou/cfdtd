#!/usr/bin/env python3
"""Draw docs/geometry/shapes/<shape>.png for every shape that has its own configuration section.

    docs/geometry/make_shape_figures.py [--fdtd ./fdtd] [--out docs/geometry/shapes] [shape ...]

Each figure is the mesh the solver builds (`fdtd mesh`, one rank, 1 nm cells) cut through the
structure, with the keys of the shape's section drawn on it as dimension lines. The dimension lines
are placed from the key values, not from the mesh, so a line that does not sit on an edge of the
drawing means the key does not mean what the figure says.

Per shape: the values used (nm, degrees for rotation) and the dimension lines, each
(view, key, (a0, b0), (a1, b1)) in nm from the centre of the structure; view "xz" is the top view
(x right, z up), "zy" the side view (z right, y up), "xy" the front view (x right, y up). `centre_cut` lists the views cut through the middle
of the structure (the default is the plane that holds the most of it). A view given
as (view, half_width) is a close-up of the centre; a line whose view is "xz:zoom" is drawn only there.
`cell` is the cell size in nm (1 by default), `n` and `ny` the cells along x, z and along y.
"""
import argparse
import os
import subprocess
import sys
import tempfile
import textwrap

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANGLES = ("rotation", "angle")  # keys in degrees
N, NY = 96, 56  # cells along x and z, along y, 1 nm cells
MATERIAL = "gold_d2cp"
STRUCTURE = 1  # the slot of the structure's material in the material maps: vacuum is 0

INI = """[grid]
nx = {n}
ny = {ny}
nz = {n}
dx = {cell}e-9
dy = {cell}e-9
dz = {cell}e-9
pml_x = 8
pml_y = 8
pml_z = 8
time_steps = 10
[spectrum]
lambda_min = 400e-9
lambda_max = 700e-9
points = 2
[source]
allow_truncated_signal = true
[geometry]
shape = {shape}
material = {material}
{geometry}
[{shape}]
{keys}
"""


def corner(cx, cz, r, sx, sz):
    """corner_radius: from the centre (cx, cz) of the corner arc to the arc, towards (sx, sz)"""
    return ("xz", "corner_radius", (cx, cz), (cx + sx * r * 0.7071, cz + sz * r * 0.7071))


def bar(p, arm_z=0.0):
    """width, length and corner radius of one rounded bar centred at z = arm_z"""
    w, l, r = p["width"], p["length"], p["corner_radius"]
    return [("xz", "width", (-w / 2, arm_z + l / 2 + 4), (w / 2, arm_z + l / 2 + 4)),
            ("xz", "length", (w / 2 + 5, arm_z - l / 2), (w / 2 + 5, arm_z + l / 2)),
            corner(-w / 2 + r, arm_z + l / 2 - r, r, -1, 1)]


def round_dipole(p):
    w, g, lt, lb, r = p["width"], p["gap"], p["top_length"], p["bottom_length"], p["corner_radius"]
    return [("xz", "width", (-w / 2, g / 2 + lt + 4), (w / 2, g / 2 + lt + 4)),
            ("xz", "top_length", (w / 2 + 5, g / 2), (w / 2 + 5, g / 2 + lt)),
            ("xz", "bottom_length", (w / 2 + 5, -g / 2), (w / 2 + 5, -g / 2 - lb)),
            ("xz", "gap", (-w / 2 - 5, -g / 2), (-w / 2 - 5, g / 2)),
            corner(-w / 2 + r, g / 2 + lt - r, r, -1, 1)]


def dipole_ends(p):
    w, g, t, b = p["width"], p["gap"], p["top_end"], p["bottom_end"]
    return [("xz", "width", (-w / 2, t + w / 2 + 4), (w / 2, t + w / 2 + 4)),
            ("xz", "gap", (-w / 2 - 5, -g / 2), (-w / 2 - 5, g / 2)),
            ("xz", "top_end", (w / 2 + 5, 0), (w / 2 + 5, t)),
            ("xz", "bottom_end", (w / 2 + 5, 0), (w / 2 + 5, -b))]


def round_hybrid(p):
    w, g, a, r = p["width"], p["gap"], p["arm_length"], p["corner_radius"]
    return [("xz", "width", (-w / 2, g / 2 + a + 4), (w / 2, g / 2 + a + 4)),
            ("xz", "arm_length", (w / 2 + 5, g / 2), (w / 2 + 5, g / 2 + a)),
            ("xz", "gap", (-w / 2 - 5, -g / 2), (-w / 2 - 5, g / 2)),
            corner(-w / 2 + r, g / 2 + a - r, r, -1, 1)]


def anapole(p):
    w = p["patch_edge"]
    r = p["corner_radius"]
    gl = p["gap_length"]
    gw = p["gap_width"]
    tip = w * (LOBE_OFFSET + LOBE_LENGTH)
    return [("xz", "patch_edge", (-w / 2, -w / 2 - 8), (w / 2, -w / 2 - 8)),
            corner(-w / 2 + r, w / 2 - r, r, -1, 1),
            ("xz", "patch_edge (lobe_offset + lobe_length)", (0, -52), (tip, -52)),
            ("xz:zoom", "gap_length", (-gl / 2, -gw / 2 - 3), (gl / 2, -gw / 2 - 3)),
            ("xz:zoom", "gap_width", (-gl / 2 - 3, -gw / 2), (-gl / 2 - 3, gw / 2))]


LOBE_LENGTH = 0.433641  # the default of [anapole] lobe_length
LOBE_OFFSET = -0.120141  # the default of [anapole] lobe_offset


def u_shape(p):
    w, l, a, r = p["width"], p["length"], p["arm_length"], p["corner_radius"]
    x0, x1 = -(a + w) / 2, (a + w) / 2
    return [("xz", "length", (x1 + 5, -l / 2), (x1 + 5, l / 2)),
            ("xz", "width", (x0, -l / 2 - 4), (x0 + w, -l / 2 - 4)),
            ("xz", "arm_length", (x0 + w, l / 2 + 4), (x1, l / 2 + 4)),
            corner(x0 + r, l / 2 - r, r, -1, 1)]


def h_shape(p):
    w, l, s, b, r = p["width"], p["length"], p["span"], p["bridge_width"], p["corner_radius"]
    return [("xz", "width", (-s / 2, l / 2 + 4), (-s / 2 + w, l / 2 + 4)),
            ("xz", "length", (s / 2 + 5, -l / 2), (s / 2 + 5, l / 2)),
            ("xz", "span", (-s / 2, -l / 2 - 4), (s / 2, -l / 2 - 4)),
            ("xz", "bridge_width", (0, -b / 2), (0, b / 2)),
            corner(s / 2 - r, l / 2 - r, r, 1, 1)]


def v_shape(p):
    w, l = p["width"], p["length"]
    c = 0.7071
    d = (w / 2 + l / 2) * c - w / (2 * c)  # x of the arm centres
    off = w / 2 + 4
    return [("xz", "length", (d - l / 2 * c + off * c, -l / 2 * c - off * c), (d + l / 2 * c + off * c, l / 2 * c - off * c)),
            ("xz", "width", (-d - w / 2 * c, -w / 2 * c), (-d + w / 2 * c, w / 2 * c))]


def l_shape(p):
    w, l, r = p["width"], p["length"], p["corner_radius"]
    return [("xz", "length", (-l / 2 - 5, -l / 2), (-l / 2 - 5, l / 2)),
            ("xz", "length", (-l / 2, -l / 2 - 4), (l / 2, -l / 2 - 4)),
            ("xz", "width", (-l / 2, l / 2 + 4), (-l / 2 + w, l / 2 + 4))]


def l_antenna(p):
    w, lx, lz = p["width"], p["length_x"], p["length_z"]
    return [("xz", "length_x", (-lx / 2, -lz / 2 - 4), (lx / 2, -lz / 2 - 4)),
            ("xz", "length_z", (-lx / 2 - 5, -lz / 2), (-lx / 2 - 5, lz / 2)),
            ("xz", "width", (-lx / 2, lz / 2 + w / 2 + 4), (-lx / 2 + w, lz / 2 + w / 2 + 4))]


def lateral_dipole(p):
    w, g, ll, lr, o = p["width"], p["gap"], p["left_length"], p["right_length"], p["overlap"]
    return [("xz", "width", (g / 2, o / 2 - lr - w / 2 - 4), (g / 2 + w, o / 2 - lr - w / 2 - 4)),
            ("xz", "gap", (-g / 2, 0), (g / 2, 0)),
            ("xz", "left_length", (-g / 2 - w - 5, -o / 2), (-g / 2 - w - 5, ll - o / 2)),
            ("xz", "right_length", (g / 2 + w + 5, o / 2), (g / 2 + w + 5, o / 2 - lr)),
            ("xz", "overlap", (-g / 2 - w - 5, -o / 2 - 16), (-g / 2 - w - 5, -o / 2 - 16 + o))]


def round_butterfly(p):
    w, g, lz, lx = p["width"], p["gap"], p["length_z"], p["length_x"]
    c = 0.7071
    return [("xz", "width", (-w / 2, w / 2 + lz + 4), (w / 2, w / 2 + lz + 4)),
            ("xz", "length_z", (w / 2 + 5, w / 2), (w / 2 + 5, w / 2 + lz)),
            ("xz", "length_x", (w / 2, -w / 2 - 5), (w / 2 + lx, -w / 2 - 5)),
            ("xz", "gap", (-g / 2 * c, -g / 2 * c), (g / 2 * c, g / 2 * c))]


def round_tripole(p):
    w, lc, lt, lb, gt, gb, r = p["width"], p["center_length"], p["top_length"], p["bottom_length"], p["top_gap"], p["bottom_gap"], p["corner_radius"]
    return [("xz", "width", (-w / 2, lc / 2 + gt + lt + 4), (w / 2, lc / 2 + gt + lt + 4)),
            ("xz", "center_length", (w / 2 + 5, -lc / 2), (w / 2 + 5, lc / 2)),
            ("xz", "top_length", (w / 2 + 5, lc / 2 + gt), (w / 2 + 5, lc / 2 + gt + lt)),
            ("xz", "bottom_length", (w / 2 + 5, -lc / 2 - gb), (w / 2 + 5, -lc / 2 - gb - lb)),
            ("xz", "top_gap", (-w / 2 - 5, lc / 2), (-w / 2 - 5, lc / 2 + gt)),
            ("xz", "bottom_gap", (-w / 2 - 5, -lc / 2), (-w / 2 - 5, -lc / 2 - gb)),
            corner(-w / 2 + r, lc / 2 + gt + lt - r, r, -1, 1)]


def x_shape(p):
    lf, wf, lr, wr, a = p["fixed_length"], p["fixed_width"], p["rotated_length"], p["rotated_width"], np.radians(p["angle"])
    u, n = np.array([np.cos(a), np.sin(a)]), np.array([-np.sin(a), np.cos(a)])
    end = (lr / 2 + 4) * u
    return [("xz", "fixed_length", (-lf / 2, -wf / 2 - 4), (lf / 2, -wf / 2 - 4)),
            ("xz", "fixed_width", (lf / 2 + 5, -wf / 2), (lf / 2 + 5, wf / 2)),
            ("xz", "rotated_length", tuple(-lr / 2 * u + (wr / 2 + 4) * n), tuple(lr / 2 * u + (wr / 2 + 4) * n)),
            ("xz", "rotated_width", tuple(end - wr / 2 * n), tuple(end + wr / 2 * n))]


SHAPES = {
    "sphere": dict(p=dict(radius=24), dims=lambda p: [("xz", "radius", (0, 0), (p["radius"], 0))]),
    "sphere_shell": dict(p=dict(radius=20), geometry="shell_material = shell\n[material.shell]\nmodel = dielectric\nn = 1.5", dims=lambda p: [("xz", "radius", (0, 0), (p["radius"], 0))]),
    "two_spheres": dict(p=dict(radius=10, separation=24), views=("xy", "xz"),
                        dims=lambda p: [("xy", "radius", (0, p["separation"] / 2), (p["radius"], p["separation"] / 2)),
                                        ("xy", "separation", (-p["radius"] - 4, -p["separation"] / 2), (-p["radius"] - 4, p["separation"] / 2))]),
    "sphere_dimer": dict(p=dict(radius=12, spacing=30),
                         dims=lambda p: [("xz", "radius", (0, p["spacing"] / 2), (p["radius"], p["spacing"] / 2)),
                                         ("xz", "spacing", (-p["radius"] - 4, -p["spacing"] / 2), (-p["radius"] - 4, p["spacing"] / 2))]),
    "sphere_trimer": dict(p=dict(radius=10, spacing=24),
                          dims=lambda p: [("xz", "radius", (0, p["spacing"]), (p["radius"], p["spacing"])),
                                          ("xz", "spacing", (-p["radius"] - 4, 0), (-p["radius"] - 4, p["spacing"]))]),
    "sphere_quadrumer": dict(p=dict(radius=8, spacing=19),
                             dims=lambda p: [("xz", "radius", (0, 1.5 * p["spacing"]), (p["radius"], 1.5 * p["spacing"])),
                                             ("xz", "spacing", (-p["radius"] - 4, 0.5 * p["spacing"]), (-p["radius"] - 4, 1.5 * p["spacing"]))]),
    "cylinder": dict(p=dict(radius=14, length=60), views=("xy", "xz"),
                     dims=lambda p: [("xy", "radius", (0, 0), (p["radius"], 0)),
                                     ("xz", "length", (p["radius"] + 5, -p["length"] / 2), (p["radius"] + 5, p["length"] / 2))]),
    "cylinder_shell": dict(p=dict(radius_core=12, core_length=64, radius_shell=20, shell_length=40), geometry="shell_material = shell\n[material.shell]\nmodel = dielectric\nn = 1.5", ny=88, views=("xz", "xy"), centre_cut=("xz",),
                           dims=lambda p: [("xz", "radius_core", (0, 0), (p["radius_core"], 0)),
                                           ("xz", "radius_shell", (0, 0), (-0.7071 * p["radius_shell"], 0.7071 * p["radius_shell"])),
                                           ("xy", "core_length", (0, -p["core_length"] / 2), (0, p["core_length"] / 2)),
                                           ("xy", "shell_length", (p["radius_shell"] + 5, -p["shell_length"] / 2), (p["radius_shell"] + 5, p["shell_length"] / 2))]),
    "ellipse": dict(p=dict(radius_x=30, radius_z=16, thickness=10, rotation=25),
                    dims=lambda p: [("xz", "radius_x", (0, 0), (p["radius_x"] * np.cos(np.radians(p["rotation"])), p["radius_x"] * np.sin(np.radians(p["rotation"])))),
                                    ("xz", "radius_z", (0, 0), (-p["radius_z"] * np.sin(np.radians(p["rotation"])), p["radius_z"] * np.cos(np.radians(p["rotation"]))))]),
    "cylinder_3d": dict(p=dict(radius=24, height=12), thickness_key="height", dims=lambda p: [("xz", "radius", (0, 0), (p["radius"], 0))]),
    "dipole": dict(p=dict(width=16, thickness=10, gap=8, top_end=32, bottom_end=22),
                   dims=lambda p: [("xz", "width", (-p["width"] / 2, p["top_end"] + 4), (p["width"] / 2, p["top_end"] + 4)),
                                   ("xz", "gap", (-p["width"] / 2 - 5, -p["gap"] / 2), (-p["width"] / 2 - 5, p["gap"] / 2)),
                                   ("xz", "top_end", (p["width"] / 2 + 5, 0), (p["width"] / 2 + 5, p["top_end"])),
                                   ("xz", "bottom_end", (p["width"] / 2 + 5, 0), (p["width"] / 2 + 5, -p["bottom_end"]))]),
    "bowtie": dict(p=dict(width=40, thickness=10, length=60, gap=8),
                   dims=lambda p: [("xz", "width", (-p["width"] / 2, p["length"] / 2 + 4), (p["width"] / 2, p["length"] / 2 + 4)),
                                   ("xz", "length", (p["width"] / 2 + 5, -p["length"] / 2), (p["width"] / 2 + 5, p["length"] / 2)),
                                   ("xz", "gap", (-p["width"] / 4 - 5, -p["gap"] / 2), (-p["width"] / 4 - 5, p["gap"] / 2))]),
    "triangle": dict(p=dict(width=40, thickness=10, length=50),
                     dims=lambda p: [("xz", "width", (-p["width"] / 2, -p["length"] / 2 - 4), (p["width"] / 2, -p["length"] / 2 - 4)),
                                     ("xz", "length", (p["width"] / 2 + 5, -p["length"] / 2), (p["width"] / 2 + 5, p["length"] / 2))]),
    "rounded_triangle": dict(p=dict(width=44, thickness=10, length=50, corner_radius=6),
                             dims=lambda p: [("xz", "width", (-p["width"] / 2, -p["length"] / 2 - 4), (p["width"] / 2, -p["length"] / 2 - 4)),
                                             ("xz", "length", (p["width"] / 2 + 5, -p["length"] / 2), (p["width"] / 2 + 5, p["length"] / 2)),
                                             corner(0, p["length"] / 2 - p["corner_radius"], p["corner_radius"], 0, 1.4142)]),
    "rounded_monopole": dict(p=dict(width=24, thickness=10, length=50, corner_radius=8), dims=lambda p: bar(p)),
    "rounded_dipole": dict(p=dict(width=20, thickness=10, length=26, gap=10, corner_radius=6),
                           dims=lambda p: bar(p, p["gap"] / 2 + p["length"] / 2) + [("xz", "gap", (-p["width"] / 2 - 5, -p["gap"] / 2), (-p["width"] / 2 - 5, p["gap"] / 2))]),
    "rounded_dipole_offset_planes": dict(p=dict(width=20, thickness=10, length=26, gap=10, corner_radius=6),
                                         dims=lambda p: bar(p, p["gap"] / 2 + p["length"] / 2) + [("xz", "gap", (-p["width"] / 2 - 5, -p["gap"] / 2), (-p["width"] / 2 - 5, p["gap"] / 2))]),
    "round_dipole": dict(p=dict(width=16, thickness=10, gap=8, top_length=28, bottom_length=20, corner_radius=4, rotation=0), dims=round_dipole),
    "one_side_round_dipole": dict(p=dict(width=14, thickness=10, gap=8, top_end=30, bottom_end=20, rotation=0), dims=dipole_ends),
    "round_hybrid": dict(p=dict(width=16, thickness=10, gap=8, arm_length=24, corner_radius=4, rotation=0), dims=round_hybrid),
    "hole": dict(p=dict(radius=16, thickness=10), dims=lambda p: [("xz", "radius", (0, 0), (p["radius"], 0))]),
    "u_shape": dict(p=dict(width=12, thickness=10, length=50, arm_length=30, corner_radius=4), dims=u_shape),
    "h_shape": dict(p=dict(width=12, thickness=10, length=50, span=44, bridge_width=10, corner_radius=4), dims=h_shape),
    "v_shape": dict(p=dict(width=12, thickness=10, length=40, corner_radius=4), dims=v_shape),
    "l_shape": dict(p=dict(width=12, thickness=10, length=40, corner_radius=4), dims=l_shape),
    "l_antenna": dict(p=dict(width=12, thickness=10, length_x=44, length_z=56), dims=l_antenna),
    "lateral_dipole": dict(p=dict(width=12, thickness=10, gap=8, left_length=44, right_length=34, overlap=14), dims=lateral_dipole),
    "round_butterfly": dict(p=dict(width=10, thickness=10, gap=8, length_z=30, length_x=22, corner_radius=3, rotation=0), dims=round_butterfly),
    "round_tripole": dict(p=dict(width=12, thickness=10, center_length=16, top_length=22, bottom_length=12, top_gap=8, bottom_gap=5, corner_radius=3, rotation=0), dims=round_tripole),
    "x_shape": dict(p=dict(thickness=10, fixed_length=64, fixed_width=12, rotated_length=44, rotated_width=8, corner_radius=3, angle=60), dims=x_shape),
    "anapole": dict(p=dict(patch_edge=200, thickness=30, corner_radius=11, gap_length=12.0432, gap_width=15.914), n=232, ny=48,
                    views=("xz", ("xz", 20), "zy"), dims=anapole),
    "two_films_spacer": dict(p=dict(thickness=10), views=("zy",), geometry="middle_thickness = 8e-9\nmiddle_medium = spacer\nlast_medium = film\n[material.spacer]\nmodel = dielectric\nn = 1.5\n[material.film]\nmodel = dielectric\nn = 2.0", dims=lambda p: []),
}

AXES = {"xz": (0, 2, 1), "zy": (2, 1, 0), "xy": (0, 1, 2)}  # horizontal, vertical, cut, axes of the [x][y][z] volume


def mesh(fdtd, shape, spec):
    keys = "\n".join(f"{k} = {v}" if k in ANGLES else f"{k} = {v}e-9" for k, v in spec["p"].items())
    text = INI.format(n=spec.get("n", N), ny=spec.get("ny", NY), cell=spec.get("cell", 1), shape=shape, material=MATERIAL, geometry=spec.get("geometry", ""), keys=keys)
    with tempfile.TemporaryDirectory() as tmp:
        ini = os.path.join(tmp, "shape.ini")
        open(ini, "w").write(text)
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        run = subprocess.run(["mpirun", "--oversubscribe", "-np", "1", fdtd, "mesh", ini, "-o", os.path.join(tmp, "out")],
                             capture_output=True, text=True, env=env, stdin=subprocess.DEVNULL)
        if run.returncode != 0:
            sys.exit(f"{shape}: fdtd mesh failed\n{run.stdout[-2000:]}\n{run.stderr[-500:]}")
        with h5py.File(os.path.join(tmp, "out", "materials", "MaterialIndex.h5")) as f:
            return f["x"][...]


def draw(shape, spec, vol, out):
    background = np.bincount(vol.ravel()).argmax()
    solid = vol != background
    if not solid.any():
        sys.exit(f"{shape}: the mesh is empty")
    idx = np.argwhere(solid)
    lo, hi = idx.min(0), idx.max(0) + 1
    # the dimension lines are placed from the centre of the domain (x, z) and of the structure (y)
    cell = spec.get("cell", 1)
    body = np.argwhere(vol == STRUCTURE)
    centre = np.array([vol.shape[0] / 2, (body[:, 1].min() + body[:, 1].max() + 1) / 2, vol.shape[2] / 2])
    views = spec.get("views", ("xz", "zy"))
    fig, axs = plt.subplots(1, len(views), figsize=(5.2 * len(views), 5.2), squeeze=False)
    dims = spec["dims"](spec["p"])
    tkey = spec.get("thickness_key", "thickness" if "thickness" in spec["p"] else None)
    for ax, view in zip(axs[0], views):
        zoom = None
        if not isinstance(view, str):
            view, zoom = view
        h, v, cut = AXES[view]
        # cut through the plane that holds the most of the structure
        at = (vol == STRUCTURE).sum(axis=(h, v)).argmax()
        if view in spec.get("centre_cut", ()):
            at = int(centre[cut])  # through the middle of the structure
        img = np.take(vol, at, axis=cut)
        img = img if h < v else img.T  # rows = horizontal axis
        ext = [-centre[h] * cell, (vol.shape[h] - centre[h]) * cell, -centre[v] * cell, (vol.shape[v] - centre[v]) * cell]
        ax.imshow((img != background).T, origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=2.2, interpolation="nearest")
        others = (img != background) & (img != STRUCTURE)
        if others.any():
            ax.imshow(np.ma.masked_where(~others.T, others.T), origin="lower", extent=ext, cmap="Blues", vmin=0, vmax=2.5, interpolation="nearest")
        for dview, key, p0, p1 in dims:
            if dview != (view + ":zoom" if zoom else view) and not (zoom is None and dview == view):
                continue
            if key == "corner_radius":
                away = (p1[0] - p0[0], p1[1] - p0[1])
                ax.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle="->", color="C3", lw=1.4, shrinkA=0, shrinkB=0))
                ax.annotate(key, p1, xytext=(18 * np.sign(away[0]) if away[0] else 0, 14 * np.sign(away[1])), textcoords="offset points", color="C3", fontsize=11,
                            ha="right" if away[0] < 0 else ("left" if away[0] > 0 else "center"), va="bottom" if away[1] > 0 else "top")
                continue
            ax.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle="<->", color="C3", lw=1.4, shrinkA=0, shrinkB=0))
            mid = ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
            vertical = abs(p1[0] - p0[0]) < abs(p1[1] - p0[1])
            below = not vertical and mid[1] < 0
            left = vertical and mid[0] < 0
            ax.annotate(key, mid, xytext=((-4 if left else 4), 0) if vertical else (0, -4 if below else 4), textcoords="offset points", color="C3", fontsize=11,
                        ha=("right" if left else "left") if vertical else "center", va="center" if vertical else ("top" if below else "bottom"))
        if tkey and view in ("zy", "xy"):
            y0, y1 = (body[:, 1].min() - centre[1]) * cell, (body[:, 1].max() + 1 - centre[1]) * cell
            x0 = min((body[:, h].max() + 1 - centre[h]) * cell + 4 * cell, ext[1] - 14 * cell)
            ax.annotate("", xy=(x0, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="<->", color="C3", lw=1.4, shrinkA=0, shrinkB=0))
            ax.annotate(tkey, (x0, (y0 + y1) / 2), xytext=(4, 0), textcoords="offset points", color="C3", fontsize=11, va="center")
        ax.set_xlabel(f"{'xyz'[h]} [nm]")
        ax.set_ylabel(f"{'xyz'[v]} [nm]")
        ax.set_title({"xz": "top view", "zy": "side view", "xy": "front view"}[view] + (", close-up" if zoom else "") + f" ({'xyz'[cut]} = {(at - centre[cut]) * cell:+.0f} nm)", fontsize=10)
        ax.set_aspect("equal")
        pad = 14 * cell
        ax.set_xlim(max(ext[0], (lo[h] - centre[h]) * cell - pad - 10 * cell), min(ext[1], (hi[h] - centre[h]) * cell + pad + 16 * cell))
        ax.set_ylim(max(ext[2], (lo[v] - centre[v]) * cell - pad), min(ext[3], (hi[v] - centre[v]) * cell + pad + 4 * cell))
        if zoom:
            ax.set_xlim(-zoom, zoom)
            ax.set_ylim(-zoom, zoom)
    values = ", ".join(f"{k} = {v:g}" + ("°" if k in ANGLES else " nm") for k, v in spec["p"].items())
    fig.suptitle("\n".join(textwrap.wrap(f"[{shape}]   {values}", 120)), fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out, shape + ".png"), dpi=110)
    plt.close(fig)


def write_page(fdtd, path):
    """docs/geometry/Shapes.md: every shape section with its keys (from `fdtd template`) and its figure"""
    text = subprocess.run([fdtd, "template"], capture_output=True, text=True, stdin=subprocess.DEVNULL).stdout
    keys, section, doc = {}, None, ""
    for line in text.split("\n"):
        if line.startswith("["):
            section = line.strip("[]")
        elif line.startswith("; "):
            doc = line[2:]
        elif "=" in line and section in SHAPES:
            keys.setdefault(section, []).append((line.split("=")[0].strip(), doc))
            doc = ""
    with open(path, "w") as f:
        f.write("# Shapes\n\nGenerated by `docs/geometry/make_shape_figures.py`; do not edit.\n\n"
                "`[geometry] shape = <name>` selects the structure and `[<name>]` holds its dimensions. Each figure is the\n"
                "mesh the solver builds for the values in its title (1 nm cells), cut through the structure, with the\n"
                "keys drawn on it. The structure is centred on `[grid] origin_x`, `origin_z`; `[geometry] shift_y` moves it\n"
                "along y. Setting a key in the section of another shape is an error.\n")
        for name in SHAPES:
            f.write(f"\n## {name}\n\n![{name}](shapes/{name}.png)\n\n| key | meaning |\n|---|---|\n")
            for key, doc in keys.get(name, []):
                f.write(f"| `{key}` | {doc} |\n")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fdtd", default=os.path.join(here, "..", "..", "fdtd"))
    ap.add_argument("--out", default=os.path.join(here, "shapes"))
    ap.add_argument("shapes", nargs="*")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for name in a.shapes or SHAPES:
        draw(name, SHAPES[name], mesh(a.fdtd, name, SHAPES[name]), a.out)
        print(name)
    write_page(a.fdtd, os.path.join(here, "Shapes.md"))
