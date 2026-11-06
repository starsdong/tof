#!/bin/bash
# sim_status3.sh [top] : progress of the ToT-test jobs (sim_job3.sh)
TOP=${1:-t5k}
cd /gpfs01/star/pwg/dongx/tof/TOFME
date
condor_q dongx -af ClusterId JobStatus | sort | uniq -c
n=$(wc -l < sim/queue_$TOP.txt)
a=$(ls output/sim/$TOP/s*/notot/simana.root 2>/dev/null | wc -l)
b=$(ls output/sim/$TOP/s*/tot/simana.root 2>/dev/null | wc -l)
c=$(ls output/sim/$TOP/s*/slow/simana.root 2>/dev/null | wc -l)
d=$(grep -l "sim_job: done" logs/sim/${TOP}_*.out 2>/dev/null | wc -l)
echo "sim $TOP: jobs $n | notot $a | tot $b | slow $c | jobs done $d"
