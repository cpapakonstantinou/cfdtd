# The material catalog

One file per material, `<name>.ini`, holding `[material.<name>]` and its poles; one file per
tensor, `chi2/<name>.ini` and `chi3/<name>.ini`. The solver reads the files a run names at start-up,
so a material is added or changed by editing a file here, or in a directory of `[materials] catalog`
or `$FDTD_MATERIALS`, which are searched first. The format and the models are in
`docs/models/Materials.md`.

The files of this directory are the parameter sets of the former built-in table, with its values;
each names the former material code in its first line.
