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
| `spectrum: <nm> amplitude ... relative ...` | the signal amplitude at each analysed wavelength; `WEAK` below 1e-2 of the strongest |

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
