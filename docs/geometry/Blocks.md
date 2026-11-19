# Blocks

A `[block.<name>]` section is an axis-aligned box of one material. Blocks build what the three
layers of `[geometry]` cannot: substrates of finite extent, buried layers, ridges, gratings.

```ini
[block.substrate]
material = glass
y_max = -50e-9          ; everything below y = -50 nm

[block.teeth]
material = silicon
x_min = 0
x_max = 500e-9
y_min = -400e-9
y_max = 0
count = 10              ; ten copies ...
pitch_x = 1e-6          ; ... 1 um apart along x
```

| Key | Default | Meaning |
|---|---|---|
| `material` | vacuum | material name ([Materials](../models/Materials.md)) |
| `x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max` | past the domain | bounds [m], measured from the structure origin (`[grid] origin_x/y/z`) |
| `count` | 1 | number of copies |
| `pitch_x`, `pitch_y`, `pitch_z` | 0 | offset between copies [m] |

- A bound that is left out reaches past the domain on that side, so `y_max` alone is a half-space.
- Blocks are drawn in file order, before the shape of `[geometry]`: a later block overwrites an
  earlier one, and the shape overwrites the blocks.
- `min < max` is required on every axis.

`fdtd mesh` shows the result (`materials/MaterialIndex_*.h5`) in a few seconds.
`examples/grating_antenna/grating_antenna.ini` builds a waveguide, a grating and an oxide cladding
from blocks.
