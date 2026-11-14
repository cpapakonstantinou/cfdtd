#!/usr/bin/env python3
"""Rewrite the [source] and [model] keys of a configuration file to the per-type sections.

    tools/migrate_config_keys.py run.ini [more.ini ...]      (files are rewritten in place)

Old [source] keys and where they go:
    kind, pulse                      -> [source] type, signal
    polarization = 90 | N            -> [plane_wave] polarization_axis = z | x, polarization_angle = N
    center_wavelength, pulse_width,
    chirp, chirp_quadratic           -> [signal.<gaussian | gaussian_alt | mode>] (pulse_width -> width)
    max_frequency                    -> [signal.<raised_cosine | raised_cosine_derivative | rectangular>]
    cw_periods, cw_settle_periods    -> [signal.cw_pulse] periods, [signal.continuous_wave] settle_periods
    wave_factor                      -> [source] aux_length_factor
    wave_derivative, wave_dispersion,
    force_aux                        -> [wave_superposition] derivative, dispersion_correction, layered_simulation
    [plane_wave.<name>]              -> [wave.<name>]
Old [model] keys:
    electrons, diffusion, fermi_velocity          -> [electrons] model, diffusion, fermi_velocity
    pressure, convection, magnetic, split_nonlinear -> [hydrodynamic]
    temperature                                   -> [temperature] model
    save_temperature, rk_stages                   -> [two_temperature] save_frames, rk_stages
    pump_probe_delay                              -> [pump_probe] delay
    update_density, shg_source, critical_point    -> dropped (read by nothing)
Other keys:
    [nonlinear] dfg_frequency, dfg_step           -> [geometry] rotation, corner_radius (what the shapes read them as)
    [geometry] height, ald_thickness, afm_size_*; [materials] stretch_z, sheet_spacing -> dropped (read by nothing)
    [lidar], [materials] damping_*, extra_layer_*  -> dropped (the lidar_* shapes are removed)
    [geometry] width, length, ... -> [<shape>] with the shape's names (every shape has a section)
        (SHAPE_KEYS below, the same table as shape_sections() in src/config/sim_config.cc)
Keys of a signal or source type that the file does not select are dropped (they had no effect), as are
circular and circular_delay; each dropped key is reported.
"""
import re
import sys

BAND = ("raised_cosine", "raised_cosine_derivative", "rectangular")
CARRIER = ("gaussian", "gaussian_alt")


def migrate(path):
    lines = open(path).read().split("\n")
    section, source = None, {}
    for line in lines:
        m = re.match(r"\s*\[([^\]]+)\]", line)
        if m:
            section = m.group(1).strip().lower()
        elif section == "source":
            m = re.match(r"\s*([A-Za-z_0-9]+)\s*=\s*([^;#]*)", line)
            if m:
                source[m.group(1).lower()] = m.group(2).strip()
    if not ({"kind", "pulse", "polarization", "center_wavelength", "pulse_width", "max_frequency"} & set(source)) \
            and not any(re.match(r"\s*\[plane_wave\.", l) for l in lines):
        return False
    signal = source.get("pulse", "raised_cosine").lower()
    mode = signal == "mode"  # the mode was a signal; it is a source type with a Gaussian signal now
    if mode:
        signal = "gaussian"
        print(f"{path}: pulse = mode -> [source] type = mode, signal = gaussian; add a [mode] section (docs/sources/ModeSource.md)")
    kind = source.get("kind", "plane_wave").lower()
    extra = {}  # section -> list of "key = value"

    def put(sect, key, value):
        extra.setdefault(sect, []).append(f"{key} = {value}")

    def drop(key, why):
        print(f"{path}: dropped [source] {key} ({why})")

    out, section = [], None
    for line in lines:
        m = re.match(r"\s*\[([^\]]+)\]", line)
        if m:
            section = m.group(1).strip().lower()
            if section.startswith("plane_wave."):
                line = line.replace("[plane_wave.", "[wave.").replace("[PLANE_WAVE.", "[wave.")
            out.append(line)
            continue
        m = re.match(r"(\s*)([A-Za-z_0-9]+)(\s*=\s*)([^;#]*)(.*)", line) if section == "source" else None
        if not m:
            out.append(line)
            continue
        key, value, tail = m.group(2).lower(), m.group(4).strip(), m.group(5)
        if key == "kind":
            out.append(f"type = {'mode' if mode else value}{'  ' + tail if tail else ''}")
        elif key == "pulse":
            out.append(f"signal = {'gaussian' if mode else value}{'  ' + tail if tail else ''}")
        elif key == "wave_factor":
            out.append(f"aux_length_factor = {value}")
        elif key == "polarization":
            if kind != "plane_wave":
                drop(key, f"kind = {kind}")
            elif int(float(value)) == 90:
                put("plane_wave", "polarization_axis", "z")
            else:
                put("plane_wave", "polarization_axis", "x")
                put("plane_wave", "polarization_angle", int(float(value)))
        elif key in ("center_wavelength", "pulse_width", "chirp", "chirp_quadratic"):
            if signal in CARRIER and not (key.startswith("chirp") and signal == "gaussian_alt"):
                put("signal." + signal, "width" if key == "pulse_width" else key, value)
            else:
                drop(key, f"pulse = {signal}")
        elif key == "max_frequency":
            if signal in BAND:
                put("signal." + signal, key, value)
            else:
                drop(key, f"pulse = {signal}")
        elif key == "cw_periods":
            put("signal.cw_pulse", "periods", value) if signal == "cw_pulse" else drop(key, f"pulse = {signal}")
        elif key == "cw_settle_periods":
            put("signal.continuous_wave", "settle_periods", value) if signal == "continuous_wave" else drop(key, f"pulse = {signal}")
        elif key in ("wave_derivative", "wave_dispersion"):
            if kind == "wave_superposition":
                put("wave_superposition", "derivative" if key == "wave_derivative" else "dispersion_correction", value)
            else:
                drop(key, f"kind = {kind}")
        elif key == "force_aux":
            if kind == "wave_superposition":
                put("wave_superposition", "layered_simulation", value)
            elif value.lower() not in ("false", "no", "0", "off"):
                drop(key, f"kind = {kind}")
        elif key in ("circular", "circular_delay"):
            drop(key, "had no effect; circular polarisation is [plane_wave] polarization_axis = rhcp | lhcp")
        else:
            out.append(line)
    # append the new sections after the [source] section
    result, section, done = [], None, False

    def flush():
        nonlocal done
        if done:
            return
        done = True
        while result and result[-1].strip() == "":
            result.pop()
        for sect, keys in extra.items():
            result.append("")
            result.append(f"[{sect}]")
            result.extend(keys)
        result.append("")

    for line in out:
        m = re.match(r"\s*\[([^\]]+)\]", line)
        if m:
            if section == "source":
                flush()
            section = m.group(1).strip().lower()
        result.append(line)
    if section == "source":
        flush()
    open(path, "w").write("\n".join(result).rstrip("\n") + "\n")
    return True


def migrate_model(path):
    lines = open(path).read().split("\n")
    keys, out, section = {}, [], None
    for line in lines:
        m = re.match(r"\s*\[([^\]]+)\]", line)
        if m:
            section = m.group(1).strip().lower()
            if section == "model":
                out.append("@MODEL@")
                continue
        if section == "model":
            m = re.match(r"\s*([A-Za-z_0-9]+)\s*=\s*([^;#]*)", line)
            if m:
                keys[m.group(1).lower()] = m.group(2).strip()
            continue
        out.append(line)
    if "@MODEL@" not in out:
        return False
    e, t = keys.get("electrons", "local").lower(), keys.get("temperature", "none").lower()
    new = ["[electrons]", f"model = {e}"]
    for k in ("diffusion", "fermi_velocity"):
        if k in keys and e != "local":
            new.append(f"{k} = {keys[k]}")
    new += ["", "[temperature]", f"model = {t}", ""]
    hyd = [f"{k} = {keys[k]}" for k in ("pressure", "convection", "magnetic", "split_nonlinear") if k in keys]
    if e == "hydrodynamic" and hyd:
        new += ["[hydrodynamic]"] + hyd + [""]
    ttm = [f"{n} = {keys[k]}" for k, n in (("save_temperature", "save_frames"), ("rk_stages", "rk_stages")) if k in keys]
    if t in ("two_temperature", "pump", "probe") and ttm:
        new += ["[two_temperature]"] + ttm + [""]
    if t in ("pump", "probe") and "pump_probe_delay" in keys:
        new += ["[pump_probe]", f"delay = {keys['pump_probe_delay']}", ""]
    for k in ("update_density", "shg_source", "critical_point"):
        if k in keys:
            print(f"{path}: dropped [model] {k} (read by nothing)")
    i = out.index("@MODEL@")
    out[i:i + 1] = new
    open(path, "w").write("\n".join(out).rstrip("\n") + "\n")
    return True


# shape -> {pool key of [geometry]: key of [<shape>]}
_BAR = {"width": "width", "thickness": "thickness", "length": "length", "arm_length": "corner_radius"}
SHAPE_KEYS = {
    "sphere": {"width": "radius"}, "sphere_shell": {"width": "radius"},
    "two_spheres": {"width": "radius", "shift_y": "separation"},
    "sphere_dimer": {"width": "radius", "length": "spacing"}, "sphere_trimer": {"width": "radius", "length": "spacing"},
    "sphere_quadrumer": {"width": "radius", "length": "spacing"},
    "cylinder": {"width": "radius"}, "ellipse": {"width": "radius_x", "length": "radius_z", "thickness": "thickness", "rotation": "rotation"}, "cylinder_3d": {"width": "radius", "length": "height"},
    "hole": {"width": "radius", "thickness": "thickness"},
    "one_side_round_dipole": {"width": "width", "thickness": "thickness", "gap": "gap", "length": "top_end", "arm_length": "bottom_end", "rotation": "rotation"},
    "round_hybrid": {"width": "width", "thickness": "thickness", "gap": "gap", "corner_radius": "corner_radius", "rotation": "rotation"},
    "dipole": {"width": "width", "thickness": "thickness", "gap": "gap", "length": "top_end", "length_2": "bottom_end"},
    "bowtie": {"width": "width", "thickness": "thickness", "length": "length", "gap": "gap"},
    "triangle": {"width": "width", "thickness": "thickness", "length": "length"},
    "rounded_triangle": _BAR, "rounded_monopole": _BAR,
    "rounded_dipole": {**_BAR, "gap": "gap"}, "rounded_dipole_offset_planes": {**_BAR, "gap": "gap"},
    "round_dipole": {"width": "width", "thickness": "thickness", "gap": "gap", "length_2": "top_length", "arm_width": "bottom_length",
                     "corner_radius": "corner_radius", "rotation": "rotation"},
    "u_shape": {**_BAR, "length_2": "arm_length"},
    "h_shape": {**_BAR, "length_2": "span", "gap": "bridge_width"},
    "v_shape": _BAR, "l_shape": _BAR,
    "two_films_spacer": {"thickness": "thickness"},
    "l_antenna": {"width": "width", "thickness": "thickness", "length": "length_x", "length_2": "length_z"},
    "lateral_dipole": {"width": "width", "thickness": "thickness", "gap": "gap", "length": "left_length", "length_2": "right_length", "arm_length": "overlap"},
    "round_butterfly": {"width": "width", "thickness": "thickness", "gap": "gap", "length_2": "length_z", "arm_width": "length_x",
                        "corner_radius": "corner_radius", "rotation": "rotation"},
    "round_tripole": {"width": "width", "thickness": "thickness", "length": "center_length", "length_2": "top_length", "arm_width": "bottom_length",
                      "arm_length": "top_gap", "gap": "bottom_gap", "corner_radius": "corner_radius", "rotation": "rotation"},
    "x_shape": {"thickness": "thickness", "length_2": "fixed_length", "width": "fixed_width", "length": "rotated_length", "gap": "rotated_width",
                "corner_radius": "corner_radius", "rotation": "angle"},
}
POOL = ("width", "length", "thickness", "length_2", "gap", "arm_length", "arm_width", "rotation", "corner_radius")
# [deg] the key is the builder's angle minus this, the round dipoles and the tripole count from +z now, the builders from +x;
# the butterfly's bars are along x and z at 0, the builder drew them along the diagonals
ROTATED = {"round_dipole": 90, "one_side_round_dipole": 90, "round_hybrid": 90, "round_tripole": 90, "round_butterfly": 45}


def migrate_shape(path):
    lines = open(path).read().split("\n")
    section, shape = None, None
    for line in lines:
        m = re.match(r"\s*\[([^\]]+)\]", line)
        if m:
            section = m.group(1).strip().lower()
        k = re.match(r"\s*shape\s*=\s*([A-Za-z_0-9]+)", line)
        if k and section == "geometry":
            shape = k.group(1).lower()
    table = SHAPE_KEYS.get(shape)
    if not table:
        return False
    out, moved, section, last, hybrid = [], [], None, None, {}
    for n, line in enumerate(lines):
        m = re.match(r"\s*\[([^\]]+)\]", line)
        if m:
            section = m.group(1).strip().lower()
        k = re.match(r"\s*([A-Za-z_0-9]+)\s*=\s*([^;#]*)", line)
        if section == "geometry" and k and not m:
            key, value = k.group(1).lower(), k.group(2).strip()
            if key in table:
                if shape in ROTATED and key == "rotation":
                    value = f"{float(value) - ROTATED[shape]:g}"
                moved.append(f"{table[key]} = {value}")
                continue
            if shape == "round_hybrid" and key in ("length_2", "arm_width", "width", "gap"):
                hybrid[key] = float(value or 0)  # arm_length = (length_2 + arm_width + width - gap) / 2
                if key in ("length_2", "arm_width"):
                    moved.append(None)
                    continue
            if key in POOL:
                if float(value or 0) != 0:
                    print(f"{path}: dropped [geometry] {key} (shape = {shape} does not read it)")
                moved.append(None)
                continue
        out.append(line)
        if section == "geometry":
            last = len(out)
    if not moved:
        return False
    keys = [x for x in moved if x]
    if shape == "round_hybrid":
        h = lambda k: hybrid.get(k, 0.0)
        keys.append(f"arm_length = {(h('length_2') + h('arm_width') + h('width') - h('gap')) / 2:g}")
    if shape in ROTATED and not any(k.startswith("rotation") for k in keys):
        keys.append(f"rotation = {-ROTATED[shape]}")  # the old default, builder angle 0
    if keys:
        while last > 0 and out[last - 1].strip() == "":
            last -= 1
        out[last:last] = ["", f"[{shape}]"] + keys
    open(path, "w").write("\n".join(out).rstrip("\n") + "\n")
    return True


def migrate_geometry(path):
    lines = open(path).read().split("\n")
    dead = {"geometry": ("height", "ald_thickness", "afm_size_x", "afm_size_z"),
            "materials": ("stretch_z", "sheet_spacing", "damping_structure", "damping_ito_perturbed",
                          "damping_ito_unperturbed", "extra_layer_fraction", "extra_layer_material")}
    moved, out, section, changed = [], [], None, False
    for line in lines:
        m = re.match(r"\s*\[([^\]]+)\]", line)
        if m:
            section = m.group(1).strip().lower()
        if section == "lidar":  # the lidar_* shapes are removed
            changed = True
            continue
        k = re.match(r"\s*([A-Za-z_0-9]+)\s*=\s*([^;#]*)", line)
        if k and not m:
            key, value = k.group(1).lower(), k.group(2).strip()
            if key in dead.get(section, ()):
                changed = True
                continue
            if section == "nonlinear" and key in ("dfg_frequency", "dfg_step"):
                changed = True
                if float(value or 0) != 0:
                    moved.append(f"{'rotation' if key == 'dfg_frequency' else 'corner_radius'} = {value}")
                continue
        out.append(line)
    if moved:
        i = next((n for n, l in enumerate(out) if re.match(r"\s*\[geometry\]", l, re.I)), None)
        if i is None:
            out += ["", "[geometry]"] + moved
        else:
            out[i + 1:i + 1] = moved
    if changed:
        open(path, "w").write("\n".join(out).rstrip("\n") + "\n")
    return changed


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        done = [migrate(p), migrate_model(p), migrate_geometry(p), migrate_shape(p)]
        print(f"{p}: {'migrated' if any(done) else 'nothing to do'}")
