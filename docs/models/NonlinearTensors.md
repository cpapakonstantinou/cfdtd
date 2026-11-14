# Nonlinear Tensors: SHG and THG

Second- and third-harmonic generation are switched on in `[nonlinear]`. The χ² or χ³ tensor of
each region is in the same file, as a named section, or is one of the built-in presets.

```ini
[nonlinear]
shg = antenna            ; none | antenna (the structure only) | bulk (the structure and the layers)
update = polarization    ; field | polarization
chi_structure = crystal  ; tensor of the structure: the [chi2.crystal] section below
chi_first = BBO          ; tensor of the first medium: a built-in preset
; chi_middle, chi_last unset: zero

[spectrum]
lambda_min = 800e-9
lambda_max = 800e-9
points = 1
harmonic_band = second   ; also analyse 2w
spacing = frequency

[chi2.crystal]           ; d_mn in m/V (chi2_ijk = 2 d_mn); components not listed are zero
d22 = 2.1e-12
d31 = 0.15e-12
```

## Keys

| Key | Meaning |
|---|---|
| `shg`, `thg` | `none`, `antenna` (the nonlinearity acts in the structure only) or `bulk` (in the structure and the layers) |
| `update` | `field`: the nonlinear term enters the E update. `polarization`: the nonlinear polarisation P_NL enters the polarisation update |
| `chi_structure` | tensor of the structure |
| `chi_first`, `chi_middle`, `chi_last` | tensors of the three layers; unset means zero |

At least one tensor has to be named. A name with neither a section nor a preset is a configuration
error, reported by `fdtd check` with the list of presets.

## Tensor sections

- `[chi2.<name>]` holds the 18 contracted coefficients `d11` … `d16`, `d21` … `d26`, `d31` … `d36`
  in m/V. The solver uses χ²_ijk = 2 d_mn.
- `[chi3.<name>]` holds the 21 components `xxxx xxyy xyyx xyxy xxzz xzzx xzxz yyyy yyxx yxxy yxyx
  yyzz yzzy yzyz zzzz zzxx zxxz zxzx zzyy zyyz zyzy` in m²/V².
- Components that are not listed are zero. An unknown component name is an error.
- Names are case-insensitive.

## Presets

| Order | Name | Content |
|---|---|---|
| χ² | `Zero` | no second-order response (centrosymmetric and isotropic media) |
| χ² | `BBO` | β-barium borate, class 3m: d11 = 1.6, d15 = d31 = 0.15, d22 = 2.1 pm/V |
| χ² | `LiNbO3` | lithium niobate, class 3m: d15 = d31 = −5.95, d33 = −34.4 pm/V |
| χ² | `Diag` | test tensor: d11 = d22 = d33 = 1 pm/V |
| χ³ | `Zero` | no third-order response |
| χ³ | `Air` | isotropic: xxxx = 1.7e-25, the other components 0.567e-25 m²/V² |
| χ³ | `Au` | gold, m3m: xxxx = 7.6e-19, the other components 2.533e-19 m²/V² |
| χ³ | `AuDiag` | gold, xxxx = yyyy = zzzz only |
| χ³ | `SiO2` | fused silica: xxxx = 2.5e-22, the other components 0.833e-22 m²/V² |
| χ³ | `SiO2Diag` | fused silica, xxxx = yyyy = zzzz only |

A section with the name of a preset replaces that preset. The presets are defined in
`src/config/chi_tensors.cc` (values after Sutherland, *Handbook of Nonlinear Optics*, 2003, and
Boyd, *Nonlinear Optics*, 3rd ed.).

## Output

`reports/NL_Efficiency_*.txt` holds the conversion efficiency; the cross-section, DFT and multipole
outputs cover the harmonic band as well. `examples/dimer_shg.ini` is a complete configuration.
