#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
date
condor_q dongx -totals | grep query
for m in 1 2 3; do for t in $(cat filelists/mudst/tags_ov.txt); do echo "mode $m $t $(grep -a "  event " logs/mudst/ov${m}_$t.log 2>/dev/null | tail -1) $(grep -a -c "StTofMEOverlayMaker::Finish" logs/mudst/ov${m}_$t.log 2>/dev/null)"; done; done
