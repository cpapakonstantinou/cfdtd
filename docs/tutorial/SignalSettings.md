# Choosing the signal and the time steps

The source signal has to end inside the run. `fdtd signal` tells you how long the signal of a
configuration is and how many time steps the run needs, before any cores are spent:

```sh
./fdtd signal run.ini
./fdtd signal -s source.signal=gaussian -s signal.gaussian.width=10e-15 run.ini   # try settings
```

It builds the signal exactly as the solver does and samples it (no MPI, well under a second).
Every line starts with a key, so `grep '^config:'` gives the lines to put in the INI file.

| Line | Meaning |
|---|---|
| `signal: dt`, `steps_configured` | the time step of the grid and `[grid] time_steps` |
| `signal: transit_steps` | steps for light to cross the longest side of the domain in vacuum |
| `signal: peak`, `last_step_on` | where the signal peaks and the last step above the threshold (`-t`, default 1e-3 of the peak) |
| `signal: steps_minimum` | the least the solver accepts |
| `signal: steps_recommended` | signal plus two transits; resonant structures ring longer than this |
| `signal: verdict` | `covered` or `TRUNCATED` |
| `config: [grid] time_steps = N` | the line to use |
| `spectrum: <nm> amplitude ... relative ...` | the signal amplitude at each analysed wavelength, DFT over the configured steps as the monitors take it; `WEAK` marks a pump point below 1e-2 of the strongest, where the normalised linear spectra are noisy; harmonic-band points are the raw DFT, not normalised, and are marked `harmonic` |

The same test runs in two more places:

- `fdtd check run.ini` (no MPI) fails when the time steps do not cover the signal;
- `fdtd simulate` makes the test first, before it allocates the fields, and stops with the
  same message.

`[source] allow_truncated_signal = true` turns the stop into a warning (short test runs).
Continuous-wave and constant-field signals have no end and are not tested.

## Harmonic runs: the noise floor at 2f and 3f

A second- or third-harmonic signal is read from the same DFT monitors as the pump, so the pump's
own spectral power at 2f and 3f is the floor the harmonic must rise above. `fdtd signal` takes the DFT
of the signal over the run's steps, as the monitors do, and reports that power relative to the
power at the fundamental:

```
harmonic: fundamental 800.0000 nm (signal carrier); [nonlinear] shg = true, thg = false
harmonic: shg 400.0000 nm floor -85.5 dB (pump power at 2 f relative to f, 40000 steps)
harmonic: thg 266.6667 nm floor -90.5 dB (pump power at 3 f relative to f, 40000 steps)
harmonic: steps 20000 shg_floor -20.1 dB thg_floor -25.1 dB
harmonic: steps 40000 shg_floor -85.5 dB thg_floor -90.5 dB (configured)
harmonic: steps 41708 shg_floor -102.2 dB thg_floor -107.7 dB (signal above the threshold covered)
harmonic: steps 80000 shg_floor < -250 dB (double-precision rounding) thg_floor < -250 dB (double-precision rounding)
```

- The fundamental is the carrier of a Gaussian signal, otherwise the centre of `[spectrum]`.
- A pulse that is cut off by the last step has a spread spectrum: the floor rises by tens of dB
  (here -20 dB at half the steps). Run until the tail has decayed well below the level you need;
  the threshold of `-t` sets where "the signal ends" for the step recommendation, the floor keeps
  dropping with the tail beyond it.
- A steady signal (continuous wave) never ends: its floor drops slowly with the run length,
  through the shrinking side lobes of the window. Pick the steps from the table.
- The floor is a property of the pump only. Whether a harmonic is visible also depends on the
  conversion efficiency, which `fdtd signal` cannot know: compare the expected harmonic level of the
  structure with the floor.

## Harmonic runs: the bands and the noise in the harmonic band

The harmonic efficiency takes its input fluence from the time signal, so the pump band does not
change the input. It does decide two things: how much of the harmonic the harmonic band (k times
every pump point) holds, and how much of the pump's own spectrum leaks into that band. With a
`[spectrum] harmonic_band`, `fdtd signal` reports both before the run:

```
./fdtd signal -s spectrum.spacing=gauss -e 1e-9 examples/cylinder_anapole.ini

band: pump 1210.0000 to 1310.0000 nm, 228.8492 to 247.7624 THz, 11 points, gauss spacing
band: pump_fwhm 31.203 THz, 165.24 nm at the fundamental (spectral intensity of the signal)
band: pump_in_band 0.5244 (share of the input fluence; ...)
band: pump_rule_error 5.77e-15 (the 11-point rule of the band against a converged integral)
band: harmonic third 403.3333 to 436.6667 nm, 686.5476 to 743.2871 THz
band: harmonic_fwhm 54.045 THz, 31.80 nm at 3 f (spectral intensity of the signal to the power 3)
band: harmonic_in_band 0.7830 (share of the harmonic energy near 3 f; ...)
band: harmonic_rule_error 4.00e-15 (the 11-point rule of the harmonic band against a converged integral)
band: verdict NARROW: 21.7 % of the harmonic falls outside the band; widen the band or lengthen the pulse
band: the band below is centred on the carrier in frequency and holds 0.999 of the harmonic
config: [spectrum] lambda_min = 1.1394e-06
config: [spectrum] lambda_max = 1.4091e-06
band: noise_efficiency 5.084e-14 (the pump's own energy in the harmonic band over the input fluence, 30000 steps; ...)
band: rounding_efficiency 4.666e-36 (double rounding of the DFT sums alone, 30000 steps)
band: margin 42.9 dB (expected efficiency 1e-09 over the larger of the two)
```

| Line | Meaning |
|---|---|
| `pump_fwhm` | width of the signal's spectral intensity; a Gaussian of field FWHM w has 0.624 / w |
| `pump_in_band` | share of the input fluence inside the pump band; matters for the linear spectra only |
| `pump_rule_error`, `harmonic_rule_error` | the solver's rule of the band (its nodes and weights) against a converged integral of the same spectrum: whether `points` is enough |
| `harmonic_fwhm`, `harmonic_in_band` | the spectrum of the signal to the power k near k f, the harmonic of an instantaneous response flat over the band, and its share inside the harmonic band; the measured efficiency is low by this factor |
| `verdict` | `covered` when the band holds the share of `-c` (default 0.999), else `NARROW` and a band centred on the carrier that holds it |
| `noise_efficiency` | the pump's own energy in the harmonic band, from the DFT over the configured steps at the solver's nodes and weights, over the input fluence |
| `rounding_efficiency` | the floor of the double-precision DFT sums alone, 2 eps^2 dt times the harmonic band width |
| `margin` | with `-e`, the expected efficiency over the larger of the two floors |

- The noise efficiency is in the units of `harmonic_efficiency*.txt` for a linear response of 1
  (a transmission of 1, or a scattering cross section equal to the normalisation area). The pump's
  content at k f reaches the monitor through the structure's linear response there, so scale it by
  that response.
- The trade: a shorter pulse or a narrower band leaves harmonic outside the band; a wider band
  holds more of the harmonic but also more of the pump's tail, and the noise efficiency rises
  (here 5.1e-14 for 1210 to 1310 nm, 1.4e-13 for the suggested 1139 to 1409 nm). A run cut before
  the signal has decayed spreads the pump's spectrum: 2.2e-8 with 26000 steps instead of 30000.
- Gauss spacing integrates these smooth spectra to rounding with a few points; uniform spacings
  use the trapezoid rule, 1.4e-3 here with 11 points in frequency.
- The harmonic estimate assumes a response flat over the band. A resonance near k f reshapes the
  harmonic: check the outermost points of the measured harmonic spectrum as well.
