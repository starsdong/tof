#!/bin/bash
# condor job: overlay/pool mixed event with production BTOF re-run, one MuDst file, one mode.
# usage: bash scripts/ov_job.sh <tag> <mode>   (mode 1 overlay+partner tracks, 2 pool, 3 pool+partner tracks)
cd /gpfs01/star/pwg/dongx/tof/TOFME
T=$1; M=$2
bash scripts/ov_test.sh -1 filelists/mudst/f_${T}.list output/mudst/ov${M}_${T}.root logs/mudst/ov${M}_${T}.log ${M}
