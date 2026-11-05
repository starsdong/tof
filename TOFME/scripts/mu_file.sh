#!/bin/bash
# condor job: method 1 (StTofMEMuMaker emulation) on one file.  usage: bash scripts/mu_file.sh <tag>
cd /gpfs01/star/pwg/dongx/tof/TOFME
T=$1
bash scripts/sl7.sh "root4star -b -q -l 'macros/runTofMEMu.C(-1,\"filelists/mudst/f_${T}.list\",\"output/mudst/mu_${T}.root\",\"geom/btof_geom_v3.txt\",5)' >&! logs/mudst/mu_${T}.log"
grep -a "StTofMEMuMaker::Finish" logs/mudst/mu_${T}.log
