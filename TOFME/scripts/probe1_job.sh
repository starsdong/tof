#!/bin/bash
# condor job: probe pass 1 (per-cell time offsets) on the first N events of one file. usage: probe1_job.sh <tag> [N]
cd /gpfs01/star/pwg/dongx/tof/TOFME
T=$1; N=${2:-4000}
bash scripts/sl7.sh "root4star -b -q -l 'macros/probeMuBTof.C(${N},\"filelists/mudst/f_${T}.list\",\"output/mudst/probe_${T}.root\",1)' >&! logs/mudst/probe_${T}_p1.log"
grep -a "probeMuBTof: done" logs/mudst/probe_${T}_p1.log
