#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
date
condor_q dongx -totals | grep query
for r in 20060065 20073019; do f=logs/mudst/me_$r.log; echo "== $r: last event line: $(grep -a "  event " $f | tail -1) | finished: $(grep -a -c "StTofMEMuMaker::Finish" $f)"; grep -a "runTofMEMu:" $f; done
