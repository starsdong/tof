#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
: > filelists/mudst/tags_ov2.txt
for m in 1 2 3; do for t in $(cat filelists/mudst/tags_ov.txt); do echo "$t, $m" >> filelists/mudst/tags_ov2.txt; done; done
wc -l filelists/mudst/tags_ov2.txt
condor_submit scripts/ov2.sub | tail -1
