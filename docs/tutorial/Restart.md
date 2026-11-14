# Stopping and restarting a run

A run that receives `SIGTERM`, `SIGINT` or `SIGUSR1`, or reaches `[run] wallclock_limit_hours`,
finishes the time step it is in, and then

1. every rank writes its state to `<output>/state/rank_<r>.bin`;
2. the outputs are written from what has been accumulated so far;
3. the run exits with `128 + signal`.

Continue it with the same configuration, output directory and rank layout:

```sh
mpirun -n 64 fdtd simulate run.ini -o /scratch/run42 --restart
```

The run is set up as usual, loads the state and continues at the next step. The result is
the one of an uninterrupted run, bit for bit. The state files are removed when the run
reaches its last step.

With SLURM, ask for the signal ahead of the time limit so the dump has time to finish:

```sh
#SBATCH --signal=B:USR1@300      # USR1 to the batch step, 300 s before the limit
```

The state holds the fields, the CPML and source-grid state, the material polarisations and
the monitors' accumulators; its size is about that of the run's field memory. The log says
`state: written after step N` or why it was not written.

Not available yet (the run still stops cleanly and writes its outputs, but leaves no
state): the hydrodynamic model, the two-temperature model, the layered-medium wave source
and the mode source.
