#!/bin/bash
# sim_status.sh [top] : progress of the simulation jobs
TOP=${1:-h1k}
cd /gpfs01/star/pwg/dongx/tof/TOFME
date
condor_q dongx -af ClusterId JobStatus | sort | uniq -c
n=$(wc -l < sim/queue_$TOP.txt)
g=$(ls output/sim/$TOP/s*/gen/hijing_*.fzd 2>/dev/null | wc -l)
t=$(ls output/sim/$TOP/s*/tfs/simana.root 2>/dev/null | wc -l)
r=$(ls output/sim/$TOP/s*/trs/simana.root 2>/dev/null | wc -l)
d=$(grep -l "sim_job: done" logs/sim/${TOP}_*.out 2>/dev/null | wc -l)
echo "sim $TOP: jobs $n | gen $g | tfs ana $t | TpcRS ana $r | jobs done $d"
