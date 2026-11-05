#!/bin/bash
# merge probe pass-1 outputs per run into output/mudst/off_<run>.root
cd /gpfs01/star/pwg/dongx/tof/TOFME
for r in $(cat filelists/mudst/runs_big.txt); do
  [ -f output/mudst/off_${r}.root ] && continue
  bash scripts/sl7.sh "hadd -f output/mudst/off_${r}.root $(ls output/mudst/probe_${r}_raw_*_pass1.root | tr '\n' ' ')" > /dev/null
done
ls output/mudst/off_*.root | wc -l
