# Nonlinear Tensors: SHG and THG

Second- and third-harmonic generation are switched on in `[nonlinear]`. The nonlinearity belongs to
the materials: a material is nonlinear when its definition names a χ² or χ³ tensor, and then every
cell of that material takes part, wherever it is drawn. [Materials](Materials.md) describes the
definitions.

```ini
[nonlinear]
shg = true               ; the chi2 tensors of the materials
update = polarization    ; field | polarization

[geometry]
material = crystal

[material.crystal]
model = dielectric
n = 1.6
chi2 = crystal           ; the [chi2.crystal] section below, or the catalog file chi2/crystal.ini

[chi2.crystal]           ; d_mn in m/V, chi2_ijk = 2 d_mn; components not listed are zero
d22 = 2.1e-12
d31 = 0.15e-12

[spectrum]
lambda_min = 800e-9
lambda_max = 800e-9
points = 1
harmonic_band = second   ; also analyse 2w
spacing = frequency
```

## Keys

| Key | Meaning |
|---|---|
| `[nonlinear] shg`, `thg` | `true`: the χ² or the χ³ tensors of the materials act; not both |
| `[nonlinear] update` | `field`: the nonlinear term enters the E update. `polarization`: the nonlinear polarisation P_NL enters the polarisation update |
| `chi2`, `chi3` of a material | the name of a `[chi2.<name>]` or `[chi3.<name>]` tensor |

A run with `shg` or `thg` needs at least one material with a tensor of that order. A material may
name both; the run loads the one of its harmonic.

A material that is linear in one place and nonlinear in another is two materials, the nonlinear one
extending the linear one:

```ini
[geometry]
material = gold_kerr     ; the structure radiates the third harmonic
last_medium = gold_d2cp  ; the gold substrate does not

[material.gold_kerr]
extends = gold_d2cp
chi3 = au
```

## Tensor sections

- `[chi2.<name>]` holds the 18 contracted coefficients `d11` … `d16`, `d21` … `d26`, `d31` … `d36`
  in m/V. The solver uses χ²_ijk = 2 d_mn.
- `[chi3.<name>]` holds the 21 components `xxxx xxyy xyyx xyxy xxzz xzzx xzxz yyyy yyxx yxxy yxyx
  yyzz yzzy yzyz zzzz zzxx zxxz zxzx zzyy zyyz zyzy` in m²/V².
- Components that are not listed are zero. An unknown component name is an error.
- A tensor is found in the run file first, then in `chi2/<name>.ini` or `chi3/<name>.ini` of the
  catalog directories.

## Tensors in the catalog

| File | Content |
|---|---|
| `chi2/zero.ini` | no second-order response, centrosymmetric and isotropic media |
| `chi2/bbo.ini` | β-barium borate, class 3m: d11 = 1.6, d15 = d31 = 0.15, d22 = 2.1 pm/V |
| `chi2/linbo3.ini` | lithium niobate, class 3m: d15 = d31 = −5.95, d33 = −34.4 pm/V |
| `chi2/diag.ini` | test tensor: d11 = d22 = d33 = 1 pm/V |
| `chi3/zero.ini` | no third-order response |
| `chi3/air.ini` | isotropic: xxxx = 1.7e-25, the other components 0.567e-25 m²/V² |
| `chi3/au.ini` | gold, m3m: xxxx = 7.6e-19, the other components 2.533e-19 m²/V² |
| `chi3/au_diag.ini` | gold, xxxx = yyyy = zzzz only |
| `chi3/sio2.ini` | fused silica: xxxx = 2.5e-22, the other components 0.833e-22 m²/V² |
| `chi3/sio2_diag.ini` | fused silica, xxxx = yyyy = zzzz only |

The values are after Sutherland, *Handbook of Nonlinear Optics*, 2003, and Boyd, *Nonlinear
Optics*, 3rd ed.

## Output

`reports/harmonic_efficiency*.txt` holds the conversion efficiency, the harmonic energy over the
input energy:

- The input fluence is n/eta0 times the integral of E_inc(t)^2 dt over the incident signal, summed
  in time over every step of the source: no band, nothing missed. n is the index of the incidence
  medium, averaged over the sampled pump spectrum when the medium is dispersive. For a Gaussian
  pulse with field FWHM w it is I0 * w * sqrt(pi / (8 ln 2)), I0 the peak intensity n E0^2 / (2 eta0).
- The harmonic exists only inside the run, so its energy comes from the DFTs: the harmonic flux at
  each point of the harmonic band, integrated over the band with the quadrature weights of the
  `[spectrum] spacing`. The harmonic band is k times every pump point, so it spans k times the pump
  band; the harmonic of a pulse is about sqrt(k) times wider than the pump, so a pump band that
  covers the pump spectrum gives a harmonic band that covers the harmonic. `spacing = gauss`
  integrates it with the fewest points.

Before 2026-11-21 the input energy came from the DFT of the source over the pump band only and both
integrals used a rectangle rule, see `docs/dev/Roadmap.md` section 7. The cross-section, DFT and multipole
outputs cover the harmonic band as well. `examples/dimer_shg.ini` is a complete configuration.
