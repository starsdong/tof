#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
bash scripts/big_status.sh
condor_q dongx -af ClusterId JobStatus | sort | uniq -c
for t in $(cat filelists/mudst/tags_all.txt); do l=logs/mudst/ov1_$t.log; [ -f $l ] && echo "$t $(grep -a "  event " $l | tail -1) $(grep -a -c StTofMEOverlayMaker::Finish $l)"; done | sort -k3 -n | tail -12
