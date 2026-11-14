#!/bin/bash
# usage: restart_test.sh name np sleep args...   (straight run vs killed + restarted run)
B=${FDTD:-./fdtd}; CMP=${CMP:-test/compare_runs.py}; name=$1; np=$2; nap=$3; shift 3
rm -rf ${OUT:-/tmp}/$name.a ${OUT:-/tmp}/$name.b
mpirun --allow-run-as-root --oversubscribe -np $np $B simulate "$@" -o ${OUT:-/tmp}/$name.a > ${OUT:-/tmp}/$name.a.log 2>&1
mpirun --allow-run-as-root --oversubscribe -np $np $B simulate "$@" -o ${OUT:-/tmp}/$name.b > ${OUT:-/tmp}/$name.b1.log 2>&1 &
sleep $nap; kill -TERM $(pgrep -f "^$B simulate" | head -1); wait $!; rc1=$?
mpirun --allow-run-as-root --oversubscribe -np $np $B simulate "$@" -o ${OUT:-/tmp}/$name.b --restart > ${OUT:-/tmp}/$name.b2.log 2>&1; rc2=$?
echo "$name: interrupted rc=$rc1 [$(grep -h '^state' ${OUT:-/tmp}/$name.b1.log | cut -c1-40)] restart rc=$rc2 [$(grep -h '^state' ${OUT:-/tmp}/$name.b2.log)] left: $(ls ${OUT:-/tmp}/$name.b/state 2>/dev/null | wc -l) | $(python3 -W ignore $CMP ${OUT:-/tmp}/$name.a ${OUT:-/tmp}/$name.b | head -4 | tr '\n' ' ')"
