#!/bin/bash
# Two-pass MuDst BTOF probe (cell offsets from matched hits, then calibrated times of all hits).
# usage: bash scripts/probe2.sh [nEvents] [list] [out]
cd /gpfs01/star/pwg/dongx/tof/TOFME
N=${1:--1}; L=${2:-filelists/mudst/test1.list}; O=${3:-output/mudst_probe.root}
for P in 1 2; do
  bash scripts/sl7.sh "root4star -b -q -l 'macros/probeMuBTof.C(${N},\"${L}\",\"${O}\",${P})' >&! logs/mudst/probe_p${P}.log"
  grep -a -E "probeMuBTof" logs/mudst/probe_p${P}.log | tail -2
done
