#!/bin/bash
# condor job: MuDst mixed event for one run.  usage: bash scripts/run_mu.sh <run>
cd /gpfs01/star/pwg/dongx/tof/TOFME
R=$1
bash scripts/sl7.sh "root4star -b -q -l 'macros/runTofMEMu.C(-1,\"filelists/mudst/run_${R}.list\",\"output/mudst/me_${R}.root\",\"geom/btof_geom_v3.txt\",5)' >&! logs/mudst/me_${R}.log"
grep -a "StTofMEMuMaker::Finish" logs/mudst/me_${R}.log
