#!/bin/bash
# per-run offset files and file lists for the MuDst mixed event (runs with completed probe pass 1)
cd /gpfs01/star/pwg/dongx/tof/TOFME
condor_q dongx -totals | grep query
: > filelists/mudst/me_runs.txt
for r in 20060065 20073019; do
  n=$(grep -a -l "probeMuBTof: done pass 1" logs/mudst/probe_${r}_raw_*_p1.log | wc -l)
  echo "run $r pass-1 done files: $n"
  hadd -f output/mudst/off_${r}.root output/mudst/probe_${r}_raw_*_pass1.root > /dev/null
  ls mudst/st_physics_${r}_raw_*.MuDst.root | sed "s|^|$PWD/|" > filelists/mudst/run_${r}.list
  echo $r >> filelists/mudst/me_runs.txt
done
ls -la output/mudst/off_*.root
