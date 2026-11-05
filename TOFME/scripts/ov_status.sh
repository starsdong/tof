#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
date
condor_q dongx -totals | grep query
for t in $(cat filelists/mudst/tags_ov.txt); do echo "$t $(grep -a "  event " logs/mudst/ov_$t.log | tail -1) $(grep -a -c "StTofMEOverlayMaker::Finish" logs/mudst/ov_$t.log)"; done
