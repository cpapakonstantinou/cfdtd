# 3D-FDTD User Guide

3D-FDTD is a parallel finite-difference time-domain solver for nano-optics. A run is described by
one INI file and started with `mpirun -n <ranks> fdtd simulate run.ini -o <output>`.

| To | Read |
|---|---|
| build the code and make a first run | [Running a Simulation](tutorial/CodeTutorial.md) |
| write a configuration file | [The Configuration File](tutorial/Controls.md) |
| choose the signal and the number of time steps | [Signal and Time Steps](tutorial/SignalSettings.md) |
| read the results | [Viewing the Outputs](tutorial/Outputs.md) |
| stop and continue a run | [Stopping and Restarting](tutorial/Restart.md) |
| choose a source | [Source Types](sources/SourceTypes.md) |
| choose a material and a material model | [Materials](models/Materials.md) |
| describe the structure | [Shapes](geometry/Shapes.md), [Blocks](geometry/Blocks.md) |
| convert an input file of the old format | [Migrating Old Input Files](dev/Migration.md) |
| change the code | [Architecture](dev/Architecture.md), [Testing](dev/Tests.md), [Style Guide](dev/StyleGuide.md) |

`fdtd template` prints every configuration key with its default and a one-line description; it is
generated from the code, so it is the reference when this guide and the code disagree.
