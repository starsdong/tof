#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
date
for t in 20073019_raw_7000013 20088015_raw_3000014 20088015_raw_6500014 20060065_raw_5500010; do f=logs/mudst/probe_${t}_p1.log; echo "== $t p1 size $(stat -c %s $f) mod $(stat -c %y $f | cut -c12-19)"; grep -a -v GARFIELD $f | tail -2 | cut -c1-150; done
condor_q dongx -run -af ClusterId ProcId RemoteHost JobCurrentStartDate | while read c p h s; do echo "$c.$p $h started $(date -d @$s +%T)"; done
