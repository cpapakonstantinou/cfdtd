# Choosing the signal and the time steps

The source signal has to end inside the run. `fdtd signal` tells you how long the signal of a
configuration is and how many time steps the run needs, before any cores are spent:

```sh
./fdtd signal run.ini
./fdtd signal -s source.signal=gaussian -s signal.gaussian.width=10e-15 run.ini
```

It builds the signal exactly as the solver does and samples it, without MPI, in well under a
second. The output is a set of tables, one key each: the first line of a table is the header, the
next lines the rows, so `grep '^band:' | cut -d' ' -f2-` is a CSV file. `config:` comes last and
holds the settings to put in the INI file.

```
./fdtd signal -s spectrum.spacing=gauss -e 1e-9 examples/cylinder_anapole.ini

signal: type,dt_s,steps,transit_steps,peak,peak_step,last_on_step,steps_minimum,steps_recommended,verdict
signal: gaussian,3.335641e-18,30000,464,2.273056e+08,18268,27240,27785,28720,covered
spectrum: lambda_nm,band,amplitude,relative,flag
spectrum: 1308.8225,pump,2.1761e-06,8.9408e-01,ok
...
spectrum: 436.2742,harmonic,4.4833e-13,1.8420e-07,raw
...
harmonic: fundamental_nm,reference
harmonic: 1260.0000,carrier
floor: steps,shg_db,thg_db,note
floor: 15000,-15.9,-21.1,
floor: 27785,-92.7,-97.4,signal_end
floor: 30000,-128.6,-135.3,configured
floor: 60000,-245.3,-249.6,
floor: 120000,-245.3,-249.6,
band: band,lambda_min_nm,lambda_max_nm,f_min_thz,f_max_thz,fwhm_thz,fwhm_nm,in_band,rule_error,verdict
band: pump,1210.0000,1310.0000,228.8492,247.7624,31.203,165.24,0.5244,5.77e-15,
band: third,403.3333,436.6667,686.5476,743.2871,54.045,31.80,0.7830,4.00e-15,narrow
noise: steps,noise_efficiency,rounding_efficiency,expected_efficiency,margin_db
noise: 30000,5.084e-14,4.666e-36,1e-09,42.9
config: section,key,value
config: grid,time_steps,30000
config: spectrum,lambda_min,1.139442e-06
config: spectrum,lambda_max,1.409088e-06
```

Warnings follow as `warning:` lines: pump points below 1e-2 of the strongest, a signal cut off by
the last step.

## The signal and the time steps

| Column of `signal:` | Meaning |
|---|---|
| `dt_s`, `steps` | the time step of the grid and `[grid] time_steps` |
| `transit_steps` | steps for light to cross the longest side of the domain in vacuum |
| `peak`, `peak_step`, `last_on_step` | the peak and the last step above the threshold, `-t`, default 1e-3 of the peak |
| `steps_minimum` | the least the solver accepts |
| `steps_recommended` | signal plus two transits; resonant structures ring longer than this |
| `verdict` | `covered`, `truncated`, or `steady` for a signal without an end |

`spectrum:` has the signal amplitude at each analysed wavelength, the DFT over the configured
steps as the monitors take it. `weak` marks a pump point below 1e-2 of the strongest, where the
normalised linear spectra are noisy; harmonic-band points are `raw`, as the solver keeps the raw
DFT there and does not normalise it.

The same test runs in two more places:

- `fdtd check run.ini`, without MPI, fails when the time steps do not cover the signal;
- `fdtd simulate` makes the test first, before it allocates the fields, and stops with the
  same message.

`[source] allow_truncated_signal = true` turns the stop into a warning, for short test runs.
Continuous-wave and constant-field signals have no end and are not tested.

## Harmonic runs: the noise floor at 2f and 3f

A second- or third-harmonic signal is read from the same DFT monitors as the pump, so the pump's
own spectral power at 2f and 3f is the floor the harmonic must rise above. `floor:` gives that power
relative to the power at the fundamental, from the DFT of the signal over a number of steps, as the
monitors take it: half, once, twice and four times the configured steps, and the end of the signal.
-250 dB stands for double-precision rounding and below.

- The fundamental is the carrier of a Gaussian signal, otherwise the centre of `[spectrum]`.
- A pulse cut off by the last step has a spread spectrum: the floor rises by tens of dB, here
  -21 dB at half the steps. Run until the tail has decayed well below the level you need; `-t` sets
  where the signal ends for the step recommendation, the floor keeps dropping with the tail beyond it.
- A steady signal never ends: its floor drops slowly with the run length, through the shrinking
  side lobes of the window. Pick the steps from the table.
- The floor is a property of the pump only. Whether a harmonic is visible also depends on the
  conversion efficiency of the structure: compare its expected harmonic level with the floor.

## Harmonic runs: the bands and the noise in the harmonic band

The harmonic efficiency takes its input fluence from the time signal, so the pump band does not
change the input. It decides two things: how much of the harmonic the harmonic band, k times every
pump point, holds, and how much of the pump's own spectrum leaks into that band. With a
`[spectrum] harmonic_band`, `band:` has a row for the pump and one for the harmonic, and `noise:`
the noise of the harmonic band.

| Column of `band:` | Meaning |
|---|---|
| `fwhm_thz`, `fwhm_nm` | width of the spectral intensity: of the signal for the pump, of the signal to the power k near k f for the harmonic; a Gaussian of field FWHM w has 0.624 / w |
| `in_band` | pump: share of the input fluence inside the band, which matters for the linear spectra only; harmonic: share of the harmonic energy inside the band, by which the measured efficiency is low |
| `rule_error` | the solver's rule of the band, its nodes and weights, against a converged integral of the same spectrum: whether `points` is enough |
| `verdict` | `covered` when the band holds the share of `-c`, default 0.999, else `narrow`, and `config:` gives a band centred on the carrier that holds it |

| Column of `noise:` | Meaning |
|---|---|
| `noise_efficiency` | the pump's own energy in the harmonic band, from the DFT over the configured steps at the solver's nodes and weights, over the input fluence |
| `rounding_efficiency` | the floor of the double-precision DFT sums alone, 2 eps^2 dt times the harmonic band width |
| `expected_efficiency`, `margin_db` | with `-e`, the expected efficiency and its margin over the larger of the two floors |

- The noise efficiency is in the units of `harmonic_efficiency*.txt` for a linear response of 1:
  a transmission of 1, or a scattering cross section equal to the normalisation area. The pump's
  content at k f reaches the monitor through the structure's linear response there, so scale it by
  that response.
- The noise adds to the harmonic as a field, not as an energy: a margin of 20 dB is an amplitude
  ratio of 10, which can move the harmonic energy by 20 %. Aim for 40 dB or more.
- The trade: a shorter pulse or a narrower band leaves harmonic outside the band; a wider band
  holds more of the harmonic but also more of the pump's tail, and the noise efficiency rises,
  here 5.1e-14 for 1210 to 1310 nm and 1.4e-13 for the suggested 1139 to 1409 nm. A run cut before
  the signal has decayed spreads the pump's spectrum: 2.2e-8 with 26000 steps instead of 30000.
- Gauss spacing integrates these smooth spectra to rounding with a few points; uniform spacings
  use the trapezoid rule, 1.4e-3 here with 11 points in frequency.
- The harmonic estimate assumes a response flat over the band. A resonance near k f reshapes the
  harmonic: check the outermost points of the measured harmonic spectrum as well.
