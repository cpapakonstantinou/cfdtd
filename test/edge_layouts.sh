#!/bin/bash
# Rank edges exactly on the CPML interfaces (base case: 52 cells per axis, interfaces at 8 and 44):
# x, y, z alone and all three, each compared with the 1-rank run.
# usage: test/edge_layouts.sh BIN TAG [extra -s options]   (run from the directory holding base.ini;
# CMP = path of test/compare_runs.py, default: next to this script)
CMP=${CMP:-$(dirname "$(readlink -f "$0")")/compare_runs.py}
B=$1; T=$2; shift 2; R=/tmp/edge_$T; mkdir -p $R
run(){ n=$1; np=$2; shift 2; rm -rf $R/$n; timeout 900 mpirun --allow-run-as-root --oversubscribe -np $np $B simulate base.ini "$@" -o $R/$n < /dev/null > $R/$n.log 2>&1 || echo "$n rc=$?"; }
run 111 1 "$@"
run x 3 -s parallel.ranks_x=3 -s "parallel.chunks_x=0 8 44" "$@"
run y 3 -s parallel.ranks_y=3 -s "parallel.chunks_y=0 8 44" "$@"
run z 3 -s parallel.ranks_z=3 -s "parallel.chunks_z=0 8 44" "$@"
run xyz 8 -s parallel.ranks_x=2 -s parallel.ranks_y=2 -s parallel.ranks_z=2 -s "parallel.chunks_x=0 44" -s "parallel.chunks_y=0 8" -s "parallel.chunks_z=0 44" "$@"
for n in x y z xyz; do echo "$T $n vs 111: $(python3 -W ignore $CMP $R/111 $R/$n 2>&1 | tail -1)"; done
