#!/bin/bash
# condor job: two-pass MuDst probe on one restored file.  usage: bash scripts/probe_job.sh <tag>
cd /gpfs01/star/pwg/dongx/tof/TOFME
T=$1
echo "/gpfs01/star/pwg/dongx/tof/TOFME/mudst/st_physics_${T}.MuDst.root" > filelists/mudst/f_${T}.list
for P in 1 2; do
  bash scripts/sl7.sh "root4star -b -q -l 'macros/probeMuBTof.C(-1,\"filelists/mudst/f_${T}.list\",\"output/mudst/probe_${T}.root\",${P})' >&! logs/mudst/probe_${T}_p${P}.log"
done
grep -a "probeMuBTof: done" logs/mudst/probe_${T}_p2.log
