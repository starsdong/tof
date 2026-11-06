#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
for f in logs/sim/h50k_*.out; do grep -q "sim_job: done" $f || { echo $f; cat $f | head -5; s=$(basename $f .out | sed s/h50k_//); ls -la --time-style=+%H:%M output/sim/h50k/s$s/gen/ | tail -4; tail -c 400 output/sim/h50k/s$s/gen/sim_$s.log; }; done
condor_q dongx -af ClusterId ProcId RemoteHost JobStartDate
