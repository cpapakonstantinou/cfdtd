#!/usr/bin/env python3
"""Rewrite the material codes of a configuration file to material names, docs/Materials.md.

    tools/migrate_materials.py run.ini [more.ini ...]      (files are rewritten in place)

What changes:
    [geometry] material, first_medium, middle_medium, last_medium, shell_material,
    [anapole] aperture_material, gap_material, [block.*] material
                                  codes -> names: 0 is vacuum, 1..99 the dielectric of [materials] n<i>,
                                  100 and up the catalog file of the former code, CODES below
    [materials] n<i>              -> [material.dielectric_<i>] model = dielectric, n = <value>
    [materials] carrier_density, carrier_density_back
                                  -> carrier_density of the ITO Drude materials, 151 .. 156
    [materials] stretch_x         -> sheet_scale of graphene, 102
    [custom_material]             -> [material.custom] of the code that read it, 302 or 351
    [nonlinear] shg, thg          -> true for antenna; bulk is gone, the layers' materials need the tensor
    [nonlinear] chi_structure     -> the structure's material becomes <material>_<tensor>, extending it with
                                     chi2 = or chi3 = <tensor>; the presets are catalog tensors, PRESETS below
    [nonlinear] chi_first, chi_middle, chi_last
                                  -> dropped: the layers were linear unless shg / thg = bulk
    [electrons] diffusion         -> diffusion of the materials that read it, 100, 101, 300, 304
Each change that is not a plain rename is reported.
"""
import re
import sys

CODES = {
    100: "gold_d2cp", 101: "silver_d2cp", 102: "graphene_d4cp", 103: "gold_silver_alloy_d2cp", 104: "gold_jc_d2cp",
    105: "gold_jc_d4cp", 106: "chromium_d2cp", 150: "gold_mcpeak_drude", 151: "ito_lidar_drude", 152: "ito_lidar_drude",
    153: "ito_10p5_drude", 154: "ito_10p5_drude", 155: "ito_2p6_drude", 156: "ito_2p6_drude", 200: "alumina_lorentz",
    201: "silica_lorentz", 202: "gaas_lorentz", 203: "lorentz_test_600thz", 250: "silicon_modified_lorentz",
    251: "gaas_modified_lorentz", 252: "silicon_modified_lorentz", 300: "gold_hydrodynamic", 301: "silver_hydrodynamic",
    304: "gold_jc_hydrodynamic", 400: "magnetic_lorentz_test",
}
ITO_FILMS = ((50e-9, "ito_film_thin"), (300e-9, "ito_film_mid"), (float("inf"), "ito_film_thick"))
PRESETS = {"zero": "", "bbo": "bbo", "linbo3": "linbo3", "diag": "diag", "air": "air", "au": "au", "audiag": "au_diag",
           "sio2": "sio2", "sio2diag": "sio2_diag"}
ROLES = {("geometry", "material"), ("geometry", "first_medium"), ("geometry", "middle_medium"), ("geometry", "last_medium"),
         ("geometry", "shell_material"), ("anapole", "aperture_material"), ("anapole", "gap_material")}
DIFFUSION_CODES = (100, 101, 300, 304)


def parse(lines):
    """section -> key -> value, the first value of each key"""
    out, section = {}, None
    for line in lines:
        m = re.match(r"\s*\[([^\]]+)\]", line)
        if m:
            section = m.group(1).strip().lower()
            out.setdefault(section, {})
            continue
        m = re.match(r"\s*([A-Za-z_0-9.]+)\s*=\s*([^;#]*)", line)
        if m and section is not None:
            out[section].setdefault(m.group(1).lower(), m.group(2).strip())
    return out


def migrate(path):
    lines = open(path).read().split("\n")
    ini = parse(lines)
    notes, added = [], {}
    geometry = ini.get("geometry", {})
    middle = float(geometry.get("middle_thickness", "0") or 0)
    materials = ini.get("materials", {})
    nonlinear = ini.get("nonlinear", {})
    custom = ini.get("custom_material", {})

    def name_of(code):
        if code == 0:
            return "vacuum"
        if code < 100:
            if "n%d" % code not in materials:
                notes.append("material %d had no [materials] n%d: its index was 0" % (code, code))
            name = "dielectric_%d" % code
            added[name] = ["model = dielectric", "n = %s" % materials.get("n%d" % code, "1")]
            return name
        if code == 350:
            return next(n for limit, n in ITO_FILMS if middle < limit)
        if code in (302, 351):
            name = "custom"
            keys = {"eps_inf": "eps_inf", "plasma_frequency": "plasma_frequency", "plasma_damping": "plasma_damping"}
            body = ["model = %s" % ("hydrodynamic" if code == 302 else "drude_lorentz")]
            body += ["%s = %s" % (k, custom[v]) for k, v in keys.items() if v in custom]
            if code == 302:
                body += ["electron_density = %s" % custom.get("electron_density", "0"),
                         "fermi_velocity = %s" % custom.get("fermi_velocity", "0")]
            added[name] = body
            poles = 2 if code == 302 else 1
            for p in range(1, poles + 1):
                suffix = "" if p == 1 else "_2"
                pole = ["amplitude = %s" % custom.get("cp_amplitude" + suffix, "0"),
                        "frequency = %s" % custom.get("cp_frequency" + suffix, "0"),
                        "damping = %s" % custom.get("cp_damping" + suffix, "0")]
                if code == 302:
                    pole.append("phase = %s" % custom.get("cp_phase" + suffix, "0"))
                added["custom.pole.%d" % p] = pole
            notes.append("code %d: [custom_material] is now [material.custom]" % code)
            return name
        if code not in CODES:
            raise SystemExit("%s: material code %d has no catalog material" % (path, code))
        base = CODES[code]
        variants = []
        if code in (151, 153, 155) and "carrier_density" in materials:
            variants.append(("carrier_density", materials["carrier_density"]))
        if code in (152, 154, 156) and "carrier_density_back" in materials:
            variants.append(("carrier_density", materials["carrier_density_back"]))
        if code == 102 and "stretch_x" in materials:
            variants.append(("sheet_scale", materials["stretch_x"]))
        diffusion = ini.get("electrons", {}).get("diffusion", "0")
        if code in DIFFUSION_CODES and float(diffusion or 0) != 0:
            variants.append(("diffusion", diffusion))
        if not variants:
            return base
        name = "%s_%d" % (base, code)
        added[name] = ["extends = %s" % base] + ["%s = %s" % kv for kv in variants]
        notes.append("code %d: %s with %s" % (code, base, ", ".join(k for k, _ in variants)))
        return name

    # the structure's tensor
    shg = nonlinear.get("shg", "none").lower()
    thg = nonlinear.get("thg", "none").lower()
    order = 2 if shg not in ("", "none", "false") else 3 if thg not in ("", "none", "false") else 0
    tensor = nonlinear.get("chi_structure", "")
    if tensor and order:
        tensor = PRESETS.get(tensor.lower(), tensor.lower())
    for mode in (shg, thg):
        if mode == "bulk":
            notes.append("shg / thg = bulk: give the layers' materials a tensor, the bulk cube is gone")
    for key in ("chi_first", "chi_middle", "chi_last"):
        t = nonlinear.get(key, "")
        if t and t.lower() != "zero":
            notes.append("[nonlinear] %s = %s dropped: the layers were linear" % (key, t))

    out, section = [], None
    skip_section = False
    for line in lines:
        m = re.match(r"\s*\[([^\]]+)\]", line)
        if m:
            section = m.group(1).strip().lower()
            skip_section = section == "custom_material"
            if not skip_section:
                out.append(line)
            continue
        if skip_section:
            continue
        m = re.match(r"(\s*)([A-Za-z_0-9.]+)(\s*=\s*)([^;#]*?)(\s*[;#].*)?$", line)
        if not m:
            out.append(line)
            continue
        key = m.group(2).lower()
        value = m.group(4).strip()
        role = (section, key) in ROLES or (section.startswith("block.") and key == "material")
        if role and re.fullmatch(r"-?\d+", value):
            code = int(value)
            if code < 0:
                continue
            name = name_of(code)
            if (section, key) == ("geometry", "material") and tensor and order:
                base = name
                name = "%s_%s" % (base, tensor)
                added[name] = ["extends = %s" % base, "chi%d = %s" % (order, tensor)]
            out.append("%s%s%s%s%s" % (m.group(1), m.group(2), m.group(3), name, m.group(5) or ""))
            continue
        if section == "materials" and (re.fullmatch(r"n\d+", key) or key in ("carrier_density", "carrier_density_back", "stretch_x")):
            continue
        if section == "nonlinear" and key in ("chi_structure", "chi_first", "chi_middle", "chi_last"):
            continue
        if section == "nonlinear" and key in ("shg", "thg"):
            out.append("%s%s%s%s%s" % (m.group(1), m.group(2), m.group(3), "false" if value.lower() in ("none", "") else "true", m.group(5) or ""))
            continue
        if section == "electrons" and key == "diffusion":
            continue
        out.append(line)
    if ("geometry", "material") not in {(s, k) for s in ini for k in ini[s]} and tensor and order:
        added["vacuum_" + tensor] = ["model = dielectric", "n = 1", "chi%d = %s" % (order, tensor)]
        notes.append("the structure was vacuum: [material.vacuum_%s] carries the tensor, set [geometry] material" % tensor)
    while out and out[-1] == "":
        out.pop()
    for name, body in added.items():
        out += ["", "[material.%s]" % name] + body
    open(path, "w").write("\n".join(out) + "\n")
    for n in notes:
        print("%s: %s" % (path, n))


if __name__ == "__main__":
    for p in sys.argv[1:]:
        migrate(p)
