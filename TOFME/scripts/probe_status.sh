#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
condor_q dongx -totals | grep query
for f in logs/mudst/probe_*_p2.log; do echo "$f: $(grep -a -h "probeMuBTof: done" $f | cut -c1-90)"; done
ls -la output/mudst/
