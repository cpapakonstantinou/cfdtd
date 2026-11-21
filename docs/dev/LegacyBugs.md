# Bugs in the original 3D-FDTD

This file lists the bugs we found in the original code while reworking it. It is meant for fixing
them in the original layout by hand. Every file and line refers to commit
`92ea34d8b9d4bd4cbb6fa9a042e28cc5aa3aa5e5` ("Update OptFunctions.cpp", 2026-03-07), the upstream
program as received. To read a file at that commit: `git show 92ea34d:09_update/updateE.cpp`.

Each entry was checked against `92ea34d` by reading the code; the "Where it was fixed in the
rework" lines give the commit of the rework that fixed it and how it was verified there. Those
commits are in the history of the rework repository, tag `backup/driver-rework-before-sync`; the
reworked code is in a different layout, so use them for the reasoning and the numbers, not as
patches.

Last updated 2026-11-21.

## Debugging the original program

Build with symbols and without optimisation, so the line numbers below hold in gdb:

```
mkdir build-debug && cd build-debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-O0 -g" ..
make -j
```

The program reads its controls from `./INIDEF` (`pphinfoini.json`, `dfts.json`, the χ files) in the
directory it runs in. One rank under gdb, from the run directory:

```
gdb -q -x bugNN.gdb ./3D-FDTD.exe
```

Several ranks, one gdb per rank, each writing its own log:

```
mpirun -np 2 sh -c 'gdb -q -batch -x bugNN.gdb ./3D-FDTD.exe > gdb.$OMPI_COMM_WORLD_RANK.log 2>&1'
```

Each "gdb commands" subsection is the content of such a `bugNN.gdb` file. Where a value has to be
chosen for your case, such as a cell inside the structure, it is set with `set $name = ...` at the
top.

## Summary

| # | Bug | Effect |
|---|---|---|
| 1 | χ table shared by dispersive and dielectric materials | the χ of a metal or other dispersive structure is overwritten by the background's: no SHG/THG from it |
| 2 | Polarisation approach of the nonlinear update | the nonlinear term is amplified by 1/(Ω dt)²; blows up in Drude-type materials at moderate fields |
| 3 | Dielectric χ with off-diagonal entries | the nonlinearity of such a dielectric is dropped without a message |
| 4 | Harmonic efficiency | input energy only over the sampled band; rectangle rule; efficiency too high |
| 5 | Drude+2CP coefficients | the pole sums run over all materials: every Drude+2CP material after the first is wrong |
| 6 | Hydrodynamic Ex update | CPML and TF/SF terms scaled by the wrong coefficient; NaN in the CPML |
| 7 | Magnetic Lorentz Ey update | averages My where Mx is meant |
| 8 | CPML α grading | integer division: α = 0 at E positions |
| 9 | Lorentz permittivity | two formulas with γ and 2γ in the code |
| 10 | Two-temperature plasma frequency, films < 50 nm | the offset reads slot 17, an empty slot |
| 11 | Plane-wave superposition, Gaussian beam, oblique injection | ranks holding CPML inject NUM_CPML cells off: results depend on the rank layout |
| 12 | Periodic ART | the reflection and transmission planes share one set of DFT arrays: T = −R on one rank in y |
| 13 | Flux boxes | the near and far face of a box share one set of DFT arrays |
| 14 | Dipole sources | a dipole on a rank boundary is injected by no rank |
| 15 | Dipole sources | the last step reads one row past the signal table |
| 16 | CW pulse signal | the carrier period is read uninitialised, in the signal and its derivative |
| 17 | Near-field point file | uninitialised values without an xz or yz DFT plane |
| 18 | 1D simulation output | the last column reads past the 1D arrays |
| 19 | Nonlinear H DFT of the XZ plane | written to the XY file name: one overwrites the other |
| 20 | Memory release | the `*_NL != NULL` checks free the linear fields a second time |

Issues seen in the rework that are probably upstream but were not checked on `92ea34d` are listed
at the end.

---

## 1. χ table shared by dispersive and dielectric materials

### Description

The second- and third-order susceptibilities are stored in `CHI2[NrDispMateD][18]` and
`CHI3[NrDispMateD][21]`. Two index spaces write into these tables:

- a dispersive material, code `c >= DispMateOffsetD` (100), writes row `c - 100`;
- a dielectric, code `c < 100`, writes row `c`.

So dispersive material 100 + s and dielectric s share row s. The rows are written in order:
antenna, first, middle, second medium. A gold antenna (code 100, row 0) in vacuum (code 0, row 0)
gets its χ from the antenna block, and the first-medium block then overwrites row 0 with the
background's χ, zero. `NL[0]`, the flag that makes the cell nonlinear, stays true. The nonlinear
update then runs in every antenna cell and adds exactly zero, in both the field and the
polarisation approach. A plasmonic antenna in air or in a dielectric of code 0 to 49 therefore has
no SHG or THG of its own; what comes out is the harmonic of the dielectrics.

The opposite direction also happens: a dielectric substrate whose code equals the slot of a
dispersive structure takes the structure's χ when its block comes later.

Physics: the harmonic generated by the metal, usually the dominant source in plasmonic SHG/THG, is
missing, while the run looks normal (nonlinear flags set, nonlinear update active). In the rework
the SHG of a Drude+2CP antenna appeared and scaled with χ after the fix; linear runs and nonlinear
dielectric structures did not change.

### Files and lines affected

| File | Lines | What |
|---|---|---|
| `06_materials/materials.cpp` | 197–226 | antenna: dispersive writes `CHI2/CHI3[ANTENNAMATERIAL - DispMateOffsetD]` (201, 205), dielectric writes `CHI2/CHI3[ANTENNAMATERIAL]` (214, 218) |
| `06_materials/materials.cpp` | 271–296 | first medium, the same two writes (275, 278 / 286, 289) |
| `06_materials/materials.cpp` | 344–369 | middle medium (348 / 359) |
| `06_materials/materials.cpp` | 417–442 | second medium (421 / 432) |
| `00_headers/globvariables.h` | 591 | `REAL CHI2[NrDispMateD][18], ...` and `CHI3` next to it |
| `00_headers/externs.h` | 629 | the matching `extern` |
| `10_models/tensorNonlinear.cpp` | 356–560 | `func_Pxyz`, `func_Pxyz_PER_COMP` read `CHI2[m]`, `CHI3[m]` with the dispersive slot |
| `10_models/isotropicNonlinear.cpp` | whole file | the dielectric exact solves read `CHI2[mx]`, `CHI3[mx]` with the dielectric code |

### gdb commands

Run a case with SHG or THG on, "Antenna Medium" = 100 (gold) with a non-zero χ file, and
"First Medium" = 0.

```
set pagination off
break main
run
watch CHI3[0][0]
commands
  silent
  printf "CHI3[0][0] = %g\n", CHI3[0][0]
  bt 2
  continue
end
watch CHI2[0][0]
commands
  silent
  printf "CHI2[0][0] = %g\n", CHI2[0][0]
  bt 2
  continue
end
break ADE_Coefficients
commands
  printf "after the material set-up: NL[0] = %d, CHI3[0][0] = %g, CHI2[0][0] = %g\n", NL[0], CHI3[0][0], CHI2[0][0]
  printf "antenna %d, first medium %d\n", pids->ANTENNAMATERIAL, pids->FIRST_MEDIUM
end
continue
```

The watchpoints show the antenna value written from `materials.cpp:201/205`, then overwritten with
0 from `materials.cpp:286/289`; at `ADE_Coefficients` `NL[0]` is 1 and the χ is 0.

### Recommended fix

Give the dielectrics their own tables, `CHI2_dl[DispMateOffsetD][18]` and
`CHI3_dl[DispMateOffsetD][21]`, as `NL_dl`, `ISO_CHI2_dl` and `ISO_CHI3_dl` already are. In the four
`else if (... < DispMateOffsetD)` blocks write `CHI2_dl[code]`/`CHI3_dl[code]` instead of
`CHI2[code]`/`CHI3[code]`, and read `CHI2_dl[mx]`/`CHI3_dl[mx]` in the dielectric exact solves.
`CHI2`/`CHI3` then hold only dispersive materials, indexed by slot.

### Files affected by the fix

`00_headers/globvariables.h`, `00_headers/externs.h` (the two new tables),
`06_materials/materials.cpp` (the four dielectric blocks), `10_models/isotropicNonlinear.cpp` (every
`CHI2[`/`CHI3[` read).

Where it was fixed in the rework: commit `7e48dea` (2026-11-06).

---

## 2. Polarisation approach of the nonlinear update

### Description

With `PNL_APPROACH = 1` ("polarization update"), `Pxyz_NL_UpdateScan` subtracts the nonlinear
polarisation P_NL = ε₀ χ E…E from the stored state of the material's auxiliary differential
equation before each E update, every step:

- Drude, Drude+CP, Lorentz slots: `Px[mmx][i][j][k] -= P_NL`;
- modified Lorentz slots: `Px_CP[mmx][0][i][j][k] -= P_NL / DT` (the pole's current).

The ADE is a second-order recursion, P^{n+1} = a₁Pⁿ + a₀Pⁿ⁻¹ + (terms in E). Changing the stored Pⁿ
without Pⁿ⁻¹ changes the polarisation and its rate together, so the subtraction acts as a force
−P_NL/dt² in d²P/dt² + γ dP/dt + Ω₀²P = …. The oscillator answers that force with

  P_eff ≈ P_NL / ((Ω₀² − Ω²) dt²)

instead of P_NL. A Drude term has no restoring force (Ω₀ = 0), so the gain is 1/(Ω dt)²: 4.0 × 10⁴
at 1.26 µm with dt = 3.34 × 10⁻¹⁸ s (2 nm cells), checked by running the recursion alone with a unit
P_NL. The explicit cubic term goes unstable once 3 χ³_eff E² reaches order one: with χ³ = 2.3 ×
10⁻¹⁸ m²/V² (ITO) that is at about 2 × 10⁶ V/m instead of about 10⁹ V/m. The blow-up is sudden, a
finite-time divergence of a cubic term, not slow exponential growth.

The damping only sets how long the Drude memory lasts (1/γ); the missing restoring term is the cause.
For a Lorentz-type pole the gain is about 1/(Ω₀ dt)²: about 470 for the UV pole of Si₃N₄ (Luke
Sellmeier) at 2 nm, too strong but usually below the blow-up. The gain grows as 1/dt², so refining
the grid makes it worse: the method does not converge.

Seen in the rework: the cylinder anapole (ITO core, Si₃N₄ shell, 10 GW/cm²) gave NaN at step
10 800, at 2 % of the pulse peak. A coarse copy reproduced it, also with a 100 times weaker pulse,
also with the ITO poles fixed, and with χ³ in the ITO alone. With χ³ in the Si₃N₄ alone, or with
`PNL_APPROACH = 0`, it stayed finite.

Both approaches take a constant χ: neither is a dispersive nonlinearity. The polarisation approach
looks frequency-dependent only through this numerical gain.

### Files and lines affected

| File | Lines | What |
|---|---|---|
| `10_models/tensorNonlinear.cpp` | 38–354 | `Pxyz_NL_UpdateScan` |
| `10_models/tensorNonlinear.cpp` | 119, 123, 218, 222, 317, 321 | `Px/Py/Pz[mm][i][j][k] -= func_Pxyz…`: Drude, Drude+CP, Lorentz |
| `10_models/tensorNonlinear.cpp` | 131, 135, 230, 234, 329, 333 | `Px/Py/Pz_CP[mm][0][i][j][k] -= func_Pxyz… / DT`: modified Lorentz |
| `09_update/updateE.cpp` | 552–559 | the call, before the E scans |
| `01_read/readInSimulation.cpp` | 965, 1063 | `PNL_APPROACH` read from the control file |

### gdb commands

A nonlinear Drude+2CP structure (for example gold, code 100) with THG on and `PNL_APPROACH = 1`.
Choose a cell inside the structure on this rank, in local indices.

```
set pagination off
set $I = 40
set $J = 40
set $K = 40
break tensorNonlinear.cpp:119 if i == $I && j == $J && k == $K
commands
  silent
  set $pnl = func_Pxyz(mmx, i, j, k, 'x')
  printf "step %d  P_NL %g  Px before %g  |Ex| %g\n", timestep, $pnl, Px[mmx][i][j][k], Ex[i][j][k]
  continue
end
run
```

P_NL is of order ε₀χE³; the stored Px drifts by many times P_NL over a few hundred steps and then
diverges with Ex. Comparing with a run at `PNL_APPROACH = 0` and the same pulse shows the field
approach staying finite.

### Recommended fix

- Do not offer the polarisation approach for slots with a Drude term (`mmx < Lorentz`): send them to
  the field approach, `Exyz_NL_UpdateScan`, or reject the combination when the input is read.
- For Lorentz-type poles, inject the nonlinear polarisation as a force scaled by the pole's own
  restoring term, Ω₀² dt² P_NL per step (`BIGOME_CRIPOI[mm][0]`), so the response at Ω ≪ Ω₀ is P_NL.
  Check the THG amplitude against the field approach at low field: they must agree.
- Until then, use `PNL_APPROACH = 0`.

### Files affected by the fix

`10_models/tensorNonlinear.cpp` (`Pxyz_NL_UpdateScan`), `09_update/updateE.cpp` (which scan runs),
`01_read/readInSimulation.cpp` (reject the polarisation approach with a Drude-type nonlinear
material).

Where it is recorded in the rework: `docs/dev/Roadmap.md` §7, "Open 2026-11-20"; not yet fixed there.

---

## 3. Dielectric χ with off-diagonal entries is ignored

### Description

For a dielectric, `materials.cpp` sets `NL_dl[code]` when any χ entry is non-zero and
`ISO_CHI2_dl`/`ISO_CHI3_dl[code]` only when every off-diagonal entry is zero. "ISO" means
diagonal-only here, not isotropic. The E scans send a nonlinear dielectric cell to the exact per-cell
solve only when it is diagonal-only; otherwise it takes the plain linear update. The field approach,
`Exyz_NL_UpdateScan`, and the polarisation approach handle dispersive slots only. A dielectric whose
χ has off-diagonal entries is therefore linear, with no message.

Physics: an isotropic χ³ needs the cross terms (χxxyy = χxyyx = χxyxy = χxxxx/3 with Kleinman
symmetry). Giving a glass, a nitride or any other dielectric its correct isotropic tensor switches
its nonlinearity off; only the diagonal-only approximation, which is not isotropic, is simulated.
Found by reading the code in the rework; not yet confirmed by a run.

### Files and lines affected

| File | Lines | What |
|---|---|---|
| `06_materials/materials.cpp` | 219–224, 291–296, 364–369, 437–442 | `ISO_CHI*_dl` only for diagonal-only tensors, `NL_dl` for any |
| `09_update/updateE.cpp` | 781, 1062, 1342 | `if (m_ == nl_code && NL_dl[m] && (ISO_CHI2_dl[m] || ISO_CHI3_dl[m]))`, else the linear update |
| `10_models/tensorNonlinear.cpp` | 637, 736, 835 | the field approach takes `m >= nl_code + DispMateOffsetD` only |
| `10_models/tensorNonlinear.cpp` | 111, 210, 309 | the polarisation approach likewise |

### gdb commands

A dielectric structure (code < 100) with THG on and an isotropic χ³ file (non-zero cross terms).

```
set pagination off
break ADE_Coefficients
commands
  printf "antenna %d: NL_dl %d ISO_CHI3_dl %d CHI3 xxxx %g xxyy %g\n", pids->ANTENNAMATERIAL, NL_dl[pids->ANTENNAMATERIAL], ISO_CHI3_dl[pids->ANTENNAMATERIAL], CHI3[pids->ANTENNAMATERIAL][0], CHI3[pids->ANTENNAMATERIAL][1]
end
break Ex_Update_THG_DL
break Exyz_NL_UpdateScan
run
```

`NL_dl` is 1 and `ISO_CHI3_dl` 0; neither `Ex_Update_THG_DL` nor a dielectric branch of
`Exyz_NL_UpdateScan` is ever reached.

### Recommended fix

Either extend the exact per-cell solve to the full tensor (a 3 × 3 Newton step per cell, the three
components coupled), or add a dielectric branch to `Exyz_NL_UpdateScan` that applies the field
approach with `CE2[m]` (the dielectric E coefficient) for codes `nl_code <= m < nl_code +
DispMateOffsetD`, with `NL_dl[m] && !ISO…_dl[m]`. At the least, stop with an error when the input
gives a dielectric off-diagonal χ.

### Files affected by the fix

`10_models/tensorNonlinear.cpp` (dielectric branch), `09_update/updateE.cpp` (the else branch at
781, 1062, 1342), or `01_read/readInSimulation.cpp` / `06_materials/materials.cpp` for the error.

Where it is recorded in the rework: `docs/dev/Roadmap.md` §7, "Open 2026-11-20".

---

## 4. Harmonic efficiency: band-limited input energy and rectangle rule

### Description

The harmonic efficiency files divide the harmonic energy by the input energy. Both are sums over
the analysed frequencies, each point times one fixed step `deltafreq`:

- the input fluence is the Parseval sum of the source DFT over the pump band only. Pump spectrum
  outside `[STOPLAMBDA ± STARTLAMBDA/2]` is not counted, so the input energy is too low and the
  efficiency too high: 5 % for a 20 fs pulse at 1.26 µm over 1210–1310 nm, a factor 7.5 for a
  broadband pulse over 400–700 nm in the rework's test;
- every point gets a full step, the edge points included: a rectangle rule over a band half a step
  wider at each end. The energies depend on the number of points (4 % between 21 and 81 points);
  the efficiency less, because the bias partly cancels.

Physics: the input energy does not need a band at all. The incident field is known in time, so the
fluence is F = (n/η₀) ∫ E_inc(t)² dt; for a Gaussian of field FWHM w, F = I₀ w √(π/(8 ln 2)). Only
the harmonic needs the DFT band, because it exists only inside the run; its band is k times every
pump point, k times as wide as the pump band, and the harmonic of a pulse is about √k times wider
than the pump, so a pump band that covers the pump covers the harmonic.

### Files and lines affected

| File | Lines | What |
|---|---|---|
| `05_sources/pulseDFT.cpp` | 178–191 | `InputIntensity`/`InputPower += 2·DT²·|S_m|²·deltafreq` over the pump band |
| `13_analysis/qCoefficient.cpp` | 344–360 | `NL_Efficiency_noPBC.txt`: `TransmNLPower += 2·DT²·Q[m]·deltafreq` (352–353) |
| `13_analysis/qCoefficient.cpp` | 522–537 | `NL_Efficiency_noPBC-ART.txt` (530–531) |
| `13_analysis/art.cpp` | 226–245 | `NL_Efficiency_PBC.txt` (238–239) |

### gdb commands

A THG or SHG run with several analysis frequencies.

```
set pagination off
break pulseDFT.cpp:192
commands
  printf "band input fluence %g J/m^2 (step %g Hz, %d points)\n", InputPower, deltafreq, NRANFREQUENCIES / 2
end
run
```

Compare with the time-domain fluence (n/η₀) Σ E(t)² dt of the source signal of the same run (sum it
in a script from the written signal file), or with I₀ w √(π/(8 ln 2)) for a Gaussian. With the pump
band narrower than the pulse, the band value is lower.

### Recommended fix

- Input: accumulate `Signal(n·DT)² · DT` over the steps of the DFT in the signal loop of
  `pulseDFT.cpp`, and set `InputPower = that · INVVACIMP · n`, with n the index of the incidence
  medium (for a dispersive medium, n averaged over the pump spectrum).
- Harmonic: weight each point with a proper quadrature rule: trapezoid (half weights at the band
  edges) for uniform points, or place the points at Gauss–Legendre nodes in frequency and use their
  weights, which integrates a smooth harmonic spectrum to many digits with a few points.

### Files affected by the fix

`05_sources/pulseDFT.cpp`, `13_analysis/qCoefficient.cpp`, `13_analysis/art.cpp`, and
`05_sources/pulseDFT.cpp` lines 30–90 for Gauss–Legendre placement of the analysis frequencies.

Where it was fixed in the rework: commit `b98e202` (2026-11-21), with `[spectrum] spacing = gauss`.

---

## 5. Drude+2CP coefficients summed over all materials

### Description

The E-update coefficients of a Drude+2CP material, `C_CRIPOI_0..3[n]`, contain the sums of
`B_CRIPOI_2/3/4[n][s]` over the material's critical-point poles. The three sums are set to zero
once before the loop over all material codes and never reset, so material n sums its own poles and
the poles of every Drude+2CP material before it.

Physics: the first Drude+2CP slot (gold, code 100) is right; silver (101) and every later Drude+2CP
material have wrong E-update coefficients, equivalent to a different permittivity. A silver sphere
changed by up to 7.6 % in the scattering and 1.2 % in the absorption cross section after the fix.

### Files and lines affected

| File | Lines | What |
|---|---|---|
| `06_materials/adeCoefficients.cpp` | 30–33 | the sums declared and zeroed before the material loop |
| `06_materials/adeCoefficients.cpp` | 98–100 | `sum_B_CRIPOI_2/3/4 += B_CRIPOI_*[n][s]` |
| `06_materials/adeCoefficients.cpp` | 102–107 | `cp9`, `cp10`, `C_CRIPOI_0..3[n]` from the sums |

### gdb commands

Any run with two Drude+2CP materials, for example antenna 101 (silver) on substrate 100 (gold).

```
set pagination off
break adeCoefficients.cpp:102
commands
  silent
  printf "slot %d: sum_B2 %g, own poles %g\n", n, sum_B_CRIPOI_2, B_CRIPOI_2[n][0] + B_CRIPOI_2[n][1]
  continue
end
run
```

For the second Drude+2CP slot the running sum differs from its own poles.

### Recommended fix

Zero `sum_B_CRIPOI_2/3/4` inside `if (n < DrudeCP)`, just before the loop over `s`.

### Files affected by the fix

`06_materials/adeCoefficients.cpp`.

Where it was fixed in the rework: commit `54d23c1` (2026-11-07).

---

## 6. Hydrodynamic Ex update: CPML and TF/SF terms scaled by C1

### Description

`Nonlin_Free_Electron_Critical_Points_Update_X` adds the CPML ψ terms and the TF/SF corrections with
`C_CRIPOI_1[mmx]`. The Y and Z versions, and every other model, use `C_CRIPOI_0`. With the
coefficients of `adeCoefficients.cpp`, C0 = DT/cp10 and C1 = 1/cp10: the Ex terms are 1/DT, about
10¹⁷, times too large.

Physics: as long as the hydrodynamic material stays away from the CPML and from the TF/SF planes the
terms are zero and nothing happens. With the material in the CPML (a substrate that runs into it) or
crossing a TF/SF plane, Ex blows up to NaN.

### Files and lines affected

| File | Lines | What |
|---|---|---|
| `10_models/hydrodynamics.cpp` | 2559 | `Nonlin_Free_Electron_Critical_Points_Update_X` |
| `10_models/hydrodynamics.cpp` | 2596, 2604, 2611, 2619 | `Ex_new ± C_CRIPOI_1[mmx] * psiEx_*` |
| `10_models/hydrodynamics.cpp` | 2627, 2632 | `Ex_new ∓ C_CRIPOI_1[mmx] * INV_DY * HzUNI[...]` (TF/SF) |
| `10_models/hydrodynamics.cpp` | 2745, 2752, 2758, 2765 | the Y version, with `C_CRIPOI_0[mmy]`, for comparison |

### gdb commands

A hydrodynamic material (codes 250–299) that reaches into the CPML or crosses a TF/SF plane.

```
set pagination off
break hydrodynamics.cpp:2596
commands
  silent
  printf "step %d cell %d %d %d: C0 %g C1 %g psi %g\n", timestep, i, j, k, C_CRIPOI_0[mmx], C_CRIPOI_1[mmx], psiEx_y_N[i][j][k]
  continue
end
run
```

C1/C0 is 1/DT; Ex in those cells reaches inf/NaN within a few hundred steps.

### Recommended fix

Replace `C_CRIPOI_1[mmx]` with `C_CRIPOI_0[mmx]` on lines 2596, 2604, 2611, 2619, 2627 and 2632.

### Files affected by the fix

`10_models/hydrodynamics.cpp`.

Where it was fixed in the rework: commit `20e148b` (2026-11-06).

---

## 7. Magnetic Lorentz Ey update averages the wrong magnetisation

### Description

`UpdateMagneticLorentz_Y` averages Mx onto the Ey position as
`0.5 * Magnetization_x[i][j][k] + 0.5 * Magnetization_y[i][j + 1][k]`. The second term is the y
component.

Physics: the magnetisation coupling of Ey mixes Mx and My; wrong whenever the magnetisation has a y
component that differs from x along y.

### Files and lines affected

`10_models/magneticLorentz.cpp`, line 74 (in `UpdateMagneticLorentz_Y`, from line 5).

### gdb commands

A magnetic Lorentz structure (codes 350–399).

```
set pagination off
break magneticLorentz.cpp:74
commands
  silent
  printf "Mx(i,j,k) %g  Mx(i,j+1,k) %g  My(i,j+1,k) used %g\n", Magnetization_x[i][j][k], Magnetization_x[i][j+1][k], Magnetization_y[i][j+1][k]
  continue
end
run
```

### Recommended fix

`Mx = 0.5 * Magnetization_x[i][j][k] + 0.5 * Magnetization_x[i][j + 1][k];`

### Files affected by the fix

`10_models/magneticLorentz.cpp`.

Where it was fixed in the rework: commit `20e148b` (2026-11-06).

---

## 8. CPML: α grading at E positions uses integer division

### Description

The complex-frequency-shift α profile at the integer (E) positions is
`MAXALPHA * pow((m) / (pids->NRCELLSCPML…), EXPONENT)` with `m` and the cell count both `int`.
The division is 0 for m < N and 1 at m = N, so α is 0 everywhere at E positions except one cell. The
half-cell (H) positions use `((m)+0.5) / N`, a real division, so they are graded as intended.

Physics: the CFS α improves the absorption of evanescent and low-frequency waves in the CPML. Here
the E and H profiles disagree. The measured effect was small, at most 0.1 dB in the reflection
(dipole corner −66.6 → −66.7 dB), fields changing at the 10⁻¹⁰ level, but the profile is not the one
intended.

### Files and lines affected

`11_boundary/cmpl.cpp`, lines 941 (x), 963 and 985 (y), 1007 (z); the H positions at 1029 and
following, for comparison.

### gdb commands

```
set pagination off
break cmpl.cpp:941
commands
  silent
  printf "m %d / N %d -> alpha %g\n", m, pids->NRCELLSCPMLX, A
  continue
end
run
```

(Put the breakpoint on line 943, after `A` is set, if the compiler places 941 before the
assignment.) Every m < N prints α = 0.

### Recommended fix

`REAL(m) / REAL(pids->NRCELLSCPMLX)` (and Y, Z) on the four lines.

### Files affected by the fix

`11_boundary/cmpl.cpp`.

Where it was fixed in the rework: commit `72a85bb` (2026-11-06).

---

## 9. Lorentz permittivity: γ in one formula, 2γ in the other

### Description

For the single-pole Lorentz materials (slots `Drude <= n < Lorentz`, codes 200–249) the code has two permittivity
formulas:

- `RefractiveIndex`: ε = ε∞ + (εs − ε∞) ωp² / (ωp² − 2iωγ − ω²);
- `ComplexPermittivity`: ε = ε∞ + (εs − ε∞) ωp² / (ωp² − iωγ − ω²).

Physics: the index serves the source normalisation in a Lorentz medium and the Q normalisation, the
permittivity other analyses; one of them has twice the damping. The ADE update of the material fixes
which convention γ means; the other formula is wrong.

### Files and lines affected

`06_materials/materials.cpp`, line 1613 (in `RefractiveIndex`, from 1585) and line 1676 (in
`ComplexPermittivity`, from 1648).

### gdb commands

A run with a Lorentz material (codes 200–249) in the incidence medium and an analysis that calls
both functions (Q coefficients or ART).

```
set pagination off
break materials.cpp:1614
commands
  silent
  printf "RefractiveIndex     slot %d omega %g: Lorentz term %g %+gi (2 gamma)\n", mmm, omega, eps._M_value
  continue
end
break materials.cpp:1677
commands
  silent
  printf "ComplexPermittivity slot %d omega %g: Lorentz term %g %+gi (gamma)\n", mmm, omega, eps._M_value
  continue
end
run
```

At the same slot and ω the two Lorentz terms differ. (With libstdc++, `print eps` also shows the
complex value if `_M_value` is not accepted.)

### Recommended fix

Derive the frequency response of the Lorentz ADE (its update coefficients in `adeCoefficients.cpp`)
and make both functions use that form. Then compute the index from `ComplexPermittivity` only, so
there is one formula.

### Files affected by the fix

`06_materials/materials.cpp`.

Where it is recorded in the rework: `docs/dev/Roadmap.md` §7, "Open 2026-11-19".

---

## 10. Two-temperature model: plasma-frequency offset from slot 17

### Description

For the ITO film model with `MIDDLE_LAYER < 50e-9`, the material is set up in slot
`DrudeLorentzShifter + 0` (= `NonLinFreeElectron` = 250), but the offset of the
temperature-dependent plasma frequency is taken from slot 17:
`PlasmaFreq_Add_Term = OME_PLASMA[17] - PlasmaFrequency(DEFAULT_TEMPERATURE)`. Slot 17 is a
Drude+CP slot (code 117), empty in this set-up, so the offset is −ωp(T₀). The other thickness
branches use `OME_PLASMA[DrudeLorentzShifter + 0]` (lines 1509, 1528); the comments `//16`, `//17`
in `definitions.h` show the slot numbers before the code ranges moved to 50-wide blocks.

Physics: the plasma frequency the two-temperature model follows is offset to near zero instead of
the film's value: wrong hot-electron dynamics for films under 50 nm.

### Files and lines affected

`06_materials/materials.cpp`, line 1491 (block 1480–1495; lines 1482–1484 also assign
`OME_PLASMA[DrudeLorentzShifter + 0]` three times, only the last counts).

### gdb commands

```
set pagination off
break materials.cpp:1491
commands
  printf "slot 17: %g, film slot %d: %g\n", OME_PLASMA[17], DrudeLorentzShifter, OME_PLASMA[DrudeLorentzShifter]
  next
  printf "PlasmaFreq_Add_Term %g\n", PlasmaFreq_Add_Term
end
run
```

### Recommended fix

`PlasmaFreq_Add_Term = OME_PLASMA[DrudeLorentzShifter + 0] - PlasmaFrequency(DEFAULT_TEMPERATURE);`
as in the other branches.

### Files affected by the fix

`06_materials/materials.cpp`.

Where it was fixed in the rework: the materials catalog commit `5de209c` (2026-11-19), where the
heated material's own plasma frequency is used.

---

## 11. Plane-wave superposition, Gaussian beam and oblique injection on CPML ranks

### Description

The superposition injectors map a local index to a global position as
`I = i + NUM_CPML_CELLS_X + startindices[0]` (and J, K). On a rank that holds part of the lower
CPML, the local index already includes the CPML shift `NRSHIFTX`, so the position is off by the
CPML thickness and that rank injects its part of the beam in the wrong place.

Physics: the fields depend on the rank layout. In the rework, 8 and 18 ranks of the same case
differed by O(1), with up to 100 % of the field leaking out of the TF/SF box; after the fix they were
bit-identical and the leakage dropped to 6 × 10⁻⁵. Layouts where no CPML rank holds an injection
plane are not affected, which hides the bug in large runs.

### Files and lines affected

`05_sources/gaussianBeam.cpp` and `05_sources/obliquePlaneWave.cpp`, 24 lines each, every
`I = i + NUM_CPML_CELLS_X + startindices[0]`, `J = j + NUM_CPML_CELLS_Y + startindices[1]`,
`K = k + NUM_CPML_CELLS_Z + startindices[2]`; the first in `gaussianBeam.cpp` are lines 50, 112,
297, 359, 455, 517, 704, 766. List them with
`git grep -n "+ NUM_CPML_CELLS_[XYZ] + startindices\[" 92ea34d -- 05_sources`.

### gdb commands

A Gaussian beam or oblique plane wave on 2 × 2 × 2 ranks; compare with 1 × 1 × 1.

```
set pagination off
break gaussianBeam.cpp:50
commands
  silent
  printf "rank %d: i %d NRSHIFTX %d start %d -> I %d (should be %d)\n", RankThis, i, NRSHIFTX, startindices[0], i + NUM_CPML_CELLS_X + startindices[0], i - NRSHIFTX + NUM_CPML_CELLS_X + startindices[0]
  continue
end
run
```

On a rank with `NRSHIFTX > 0` the two positions differ.

### Recommended fix

Subtract the local shift on every such line: `I = i - NRSHIFTX + NUM_CPML_CELLS_X + startindices[0]`,
and likewise `- NRSHIFTY`, `- NRSHIFTZ`.

### Files affected by the fix

`05_sources/gaussianBeam.cpp`, `05_sources/obliquePlaneWave.cpp`.

Where it was fixed in the rework: commit `f25970b` (2026-11-03). The layered-medium injector had a
separate layout dependence at the TF/SF box corners, also in the original; see the end.

---

## 12. Periodic ART: reflection and transmission planes share one set of arrays

### Description

The ART planes of periodic runs, near y (`JART1`, reflection) and far y (`JART2`, transmission),
both accumulate into the same DFT arrays, `Re/ImTFEx/Ez/Hx/Hz_XZ_ART`. When one rank holds both
planes, which is every run with one rank along y, `Reflectance_XZ` and `Transmittance_XZ` integrate
the sum of the two planes.

Physics: T = −R in every ART file and A = 1 exactly. With the planes on different ranks the results
are right, so the bug depends on the layout.

### Files and lines affected

| File | Lines | What |
|---|---|---|
| `08_monitors/updateDFT.cpp` | 252–253, 261–262 | both planes into `*_XZ_ART` |
| `13_analysis/art.cpp` | 35, 57 | `Reflectance_XZ` and `Transmittance_XZ` on the same arrays |
| `13_analysis/art.cpp` | 285, 356 | the two integrals |
| `02_initialization/setupSimulation.cpp` | 2794 and following | one set of `*_XZ_ART` arrays allocated |

### gdb commands

A periodic run with ART on 1 rank along y.

```
set pagination off
break updateDFT.cpp:252
commands
  silent
  printf "step %d: near plane %d into _ART\n", timestep, JART1
  continue
end
break updateDFT.cpp:261
commands
  silent
  printf "step %d: far plane %d into the same _ART\n", timestep, JART2
  continue
end
run
```

Both fire on the same rank; the ART output has T = −R.

### Recommended fix

Allocate a second set, `*_XZ_ART_FAR`, accumulate `JART2` into it, and pass the near set to
`Reflectance_XZ` and the far set to `Transmittance_XZ`. The far-field point monitors then read the
far set (forward) and the near set (backward).

### Files affected by the fix

`08_monitors/updateDFT.cpp`, `13_analysis/art.cpp`, `02_initialization/setupSimulation.cpp`
(allocation), `03_memory/memoryRelease.cpp`, `08_monitors/normalizeDFT.cpp` (lines 61–64),
`08_monitors/updateFieldMonitors.cpp` (far-field points), `00_headers/globvariables.h`,
`00_headers/externs.h`.

Where it was fixed in the rework: commit `c464967` (2026-11-07).

---

## 13. Flux boxes: near and far face share one set of arrays

### Description

The scattering and total-field boxes have a near and a far face per axis (`QSCANEAX`/`QSCAFARX`,
...). Both faces of an axis accumulate into the same DFT arrays (`*_YZ_SCA`, `*_XZ_SCA`,
`*_XY_SCA`, and `*_TOT` for the total box). On a rank that holds both faces of a box along an axis,
which is every 1-rank run and every run with one rank along that axis, the two faces add up.

Physics: the flux through the box is wrong; in the rework a 1-rank run gave 3.8 times the
transmission. With the faces on different ranks it is right. A related problem in the same code: a
rank that holds a face plane but does not overlap the face transversally kept the face flag and
injected or sampled at local cell 0.

### Files and lines affected

`08_monitors/updateDFT.cpp`, lines 38–185: the pairs `QSCANEAX`/`QSCAFARX` (38, 51),
`QSCANEAY`/`QSCAFARY` (63, 75), `QSCANEAZ`/`QSCAFARZ` (87, 100), `QTOTNEA*`/`QTOTFAR*` (113–185),
each pair writing the same arrays.

### gdb commands

Any case with flux boxes on 1 rank.

```
set pagination off
break updateDFT.cpp:63
commands
  silent
  printf "step %d: near y face %d into _XZ_SCA\n", timestep, JQSCANEA
  continue
end
break updateDFT.cpp:75
commands
  silent
  printf "step %d: far y face %d into the same _XZ_SCA\n", timestep, JQSCAFAR
  continue
end
run
```

### Recommended fix

Give every far face its own arrays (`*_FAR`), integrate each face from its own set, and clear a
rank's face flag when it does not overlap the face transversally.

### Files affected by the fix

`08_monitors/updateDFT.cpp`, `13_analysis/qCoefficient.cpp`, `02_initialization/setupSimulation.cpp`
(allocation and face flags), `03_memory/memoryRelease.cpp`, `08_monitors/normalizeDFT.cpp`,
`00_headers/globvariables.h`, `00_headers/externs.h`.

Where it was fixed in the rework: commit `54e00f6` (2026-11-03), with the flux-box ranges and the
normalisation area of each face.

---

## 14. Dipole on a rank boundary is injected by no rank

### Description

`SetupDipoleSources` keeps a dipole on a rank when
`X_dipole < startindices[0] + localsize[0] - 1` (and Y, Z). A dipole whose cell lies on the plane
two ranks share fails the test on the lower rank (it is its last cell) and lies before the start of
the upper rank's own range on the other, so no rank injects it.

Physics: a 2-rank run can have no source at all where the 1-rank run has one.

### Files and lines affected

`02_initialization/setupSources.cpp`, line 487 (in `SetupDipoleSources`, from 438).

### gdb commands

One dipole on a rank boundary, for example x = N/2 on 2 ranks along x.

```
set pagination off
break setupSources.cpp:487
commands
  silent
  printf "rank %d: dipole %d at %d %d %d, start %d size %d\n", RankThis, i, X_dipole, Y_dipole, Z_dipole, startindices[0], localsize[0]
  continue
end
break setupSources.cpp:512
commands
  printf "rank %d keeps %d dipoles\n", RankThis, NUM_DIPOLE_SOURCES_THIS
  continue
end
run
```

### Recommended fix

Keep the dipole on every rank whose cell range contains it, `< startindices[a] + localsize[a]`, so
both ranks inject it into their copy of the shared cell, as the field update of a shared cell does.

### Files affected by the fix

`02_initialization/setupSources.cpp`.

Where it was fixed in the rework: commit `15db606` (2026-11-06).

---

## 15. Dipole: the last step reads past the signal table

### Description

`DipoleSource` has `NRMAXTIMESTEPS` rows (steps 0 to N − 1), but the time loop runs
`timestep = 1 … NRMAXTIMESTEPS` and reads `DipoleSource[timestep]`. The last step reads row N, one
past the table.

Physics: the last step injects whatever is in memory after the table (0 on the rework's machine,
garbage in general). Undefined behaviour.

### Files and lines affected

| File | Lines | What |
|---|---|---|
| `02_initialization/setupSources.cpp` | 512–529 | `MemAllocat3DArray(DipoleSource, NRMAXTIMESTEPS, …)`, filled for j < NRMAXTIMESTEPS |
| `main.cpp` | 1716 | `for (timestep = 1; timestep <= NRMAXTIMESTEPS; timestep++)` |
| `main.cpp` | 1040–1042, 1219–1221 | `DipoleSource[timestep][i][…]` |

### gdb commands

```
set pagination off
break main.cpp:1219 if timestep == NRMAXTIMESTEPS
commands
  printf "step %d reads row %d of a table with %d rows: %g\n", timestep, timestep, NRMAXTIMESTEPS, DipoleSource[timestep][i][0]
end
run
```

Or run under valgrind: `valgrind --track-origins=yes ./3D-FDTD.exe` reports an invalid read there.

### Recommended fix

Allocate and fill `NRMAXTIMESTEPS + 1` rows (`j <= NRMAXTIMESTEPS`).

### Files affected by the fix

`02_initialization/setupSources.cpp`, and its release in `03_memory/memoryRelease.cpp` with the new
size.

Where it was fixed in the rework: commit `345d2de` (2026-11-07).

---

## 16. CW pulse: the carrier period is read uninitialised

### Description

In `Signal` and `dSignal` the CW pulse branch computes `k_period = PI2 / period` and the carrier
delay `0.5 * num_periods_cw * period`, but `period` is a local variable of the function that this
branch never sets (only the branches of `Envelopes` set their own copy).

Physics: the carrier frequency and its delay of a `cw_pulse` signal are whatever is on the stack:
the SHG/THG CW-pulse runs are undefined. The intended period is 1 / fundamental frequency, which the
envelope and the step count already assume.

### Files and lines affected

`05_sources/signal.cpp`: `Signal` declares `period` at 156 and uses it at 171–175; `dSignal` declares
it at 432 and uses it at 450–455.

### gdb commands

A run with the CW pulse signal.

```
set pagination off
break signal.cpp:173
commands
  silent
  printf "t %g: period %g (uninitialised), expected %g\n", t, period, 1.0 / StartFrequency
  continue
end
run
```

(`StartFrequency` or whichever variable holds the fundamental frequency in the run.)

### Recommended fix

Set `period = 1 / fundamental frequency` at the top of both branches, from the same value the CW
envelope and the step count use.

### Files affected by the fix

`05_sources/signal.cpp`.

Where it was fixed in the rework: commits `34a7561` (signal, 2026-11-12) and `04e3ffc` (derivative,
2026-11-13).

---

## 17. Near-field point file written from uninitialised values

### Description

The near-field point writer sets `Ex_mod`, `Ez_mod` and the phases only inside `if (DFT_YZ)` or
`if (DFT_XZ)`. With the animation planes on but neither DFT plane, it still opens
`NEAR_FIELD_<pol>.txt` and prints the uninitialised values.

Physics: a file of garbage that looks like data.

### Files and lines affected

`08_monitors/updateFieldMonitors.cpp`, lines 128–162 (`if (near_field_point)` at 130, `if (DFT_YZ)`
at 140, `if (DFT_XZ)` at 147, the print at 160).

### gdb commands

```
set pagination off
break updateFieldMonitors.cpp:160
commands
  silent
  printf "DFT_YZ %d DFT_XZ %d -> Ex_mod %g Ez_mod %g\n", DFT_YZ, DFT_XZ, Ex_mod, Ez_mod
  continue
end
run
```

### Recommended fix

Write the file only when `DFT_YZ || DFT_XZ`, or initialise the values and say in the file that no
plane exists.

### Files affected by the fix

`08_monitors/updateFieldMonitors.cpp`.

Where it was fixed in the rework: commit `2e353c6` (2026-11-07).

---

## 18. 1D simulation output reads past the 1D arrays

### Description

The writer of `1D_SIMULATION/*.txt` loops `j = -gridShift … gridLength - gridShift - 1` and
interpolates the staggered fields at each point; at the last point the interpolation reads index
`gridLength`, one past the 1D arrays.

Physics: the last column of every 1D file is uninitialised memory and changes from run to run (Ey
of a z-polarised wave was 4 × 10⁻⁵ in one run and 1.4 in the next, 0 expected). The 3D injection
never reaches that point, so the fields are not affected.

### Files and lines affected

`05_sources/oneDimensionalSimulator.cpp`, line 67 (loop, block 55–90); the array length is set in
`05_sources/planeWaveOneDimensionalSimulator.cpp`, line 71.

### gdb commands

Run under valgrind, which reports the invalid read:

```
valgrind --track-origins=yes ./3D-FDTD.exe 2>&1 | grep -A8 "Invalid read" | head -40
```

or in gdb:

```
set pagination off
break oneDimensionalSimulator.cpp:70 if j == gridLength - gridShift - 1
commands
  printf "last point j %d -> index %d of %d\n", j, j + gridShift, gridLength
  step
end
run
```

### Recommended fix

Stop one point earlier: `j < gridLength - gridShift - 1`.

### Files affected by the fix

`05_sources/oneDimensionalSimulator.cpp`.

Where it was fixed in the rework: commit `979278c` (2026-11-07).

---

## 19. Nonlinear H DFT of the XZ plane written to the XY file name

### Description

In the XZ block of the DFT writer, the nonlinear magnetic file is named
`"/Nonlinear_DFTm_XY_" + …`, the same name as the XY block's file.

Physics: with both planes on, one file overwrites the other; the XZ nonlinear H is lost or read as
XY.

### Files and lines affected

`14_write/writeDFT.cpp`, line 336 (XZ block, which names its E file `Nonlinear_DFT_XZ_` at 312); the
XY block's own file at 205.

### gdb commands

```
set pagination off
break writeDFT.cpp:336
commands
  printf "XZ block writes %s\n", strFileH.c_str()
end
run
```

### Recommended fix

`"/Nonlinear_DFTm_XZ_"` on line 336.

### Files affected by the fix

`14_write/writeDFT.cpp`.

Where it was fixed in the rework: commit `172d30f` (2026-11-03).

---

## 20. Memory release: the `*_NL` checks free the linear fields again

### Description

`MemoryAndDataRelease` frees `Ex … Hz` (lines 25–33), then:

```
if (Ex_NL != NULL)
    MemRelease3DArray(Ex, …);
```

and the same for the other five components: the check is on the nonlinear field but the release is
of the linear one. `MemRelease3DArray` deletes the blocks but does not set the pointer to NULL, so
the `!= NULL` checks cannot catch a second release.

Physics: none on the fields, but when the nonlinear fields exist (hydrodynamic split runs) `Ex`…`Hz`
are freed twice (heap corruption, a crash at the end of the run or in a later sweep) and `Ex_NL`…
`Hz_NL` are never freed.

### Files and lines affected

| File | Lines | What |
|---|---|---|
| `03_memory/memoryRelease.cpp` | 34–45 | `if (E*_NL / H*_NL != NULL) MemRelease3DArray(E* / H*, …)` |
| `03_memory/allocateMemory.cpp` | 121–133 | `MemRelease3DArray` does not reset the pointer |

### gdb commands

A hydrodynamic split run (the nonlinear fields allocated).

```
set pagination off
break MemoryAndDataRelease
run
break MemRelease3DArray<double>
commands
  silent
  printf "release %p\n", p
  bt 1
  continue
end
continue
```

The same address is released twice. With glibc, `MALLOC_CHECK_=3 ./3D-FDTD.exe` aborts at the
second release.

### Recommended fix

Release `Ex_NL`…`Hz_NL` in those six lines, and set `p = NULL` at the end of `MemRelease3DArray` (and
the other release helpers).

### Files affected by the fix

`03_memory/memoryRelease.cpp`, `03_memory/allocateMemory.cpp`.

Where it was fixed in the rework: the arrays became owning smart pointers (commit `77844aa`,
2026-11-03), and the release of the state was rewritten on 2026-11-13 (`docs/dev/Roadmap.md` §7,
"Fixed 2026-11-13").

---

## Probably upstream, not checked on 92ea34d

These were seen in the rework and described there as pre-existing or also upstream, but were not
traced to lines of `92ea34d`:

- Two-temperature model: `Update_TTM` crashes (`docs/dev/Roadmap.md` §7).
- Magnetic Lorentz as a background medium (code 400 in the first/middle/last medium) segfaults; as
  a structure it runs. Its update also has no CPML ψ or TF/SF terms, so it is only usable inside the
  TF/SF box.
- Hydrodynamic runs with a single analysis wavelength write empty `Qcoeff` tables.
- The layered-medium (oblique) injector depends on the rank layout at the TF/SF box corners.
- The hydrodynamic LWA multipole writer had no code filling its 60 time-series arrays, and its only
  call needed `e_dens_ft`, which was always false.
- No check for non-finite fields: a run that has gone to NaN keeps running to the last step.
