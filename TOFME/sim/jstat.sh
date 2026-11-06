#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
condor_q dongx -af ClusterId JobStatus
tail -2 logs/sim/jtest_2001.out
f=output/sim/jtest/s2001/trs/reco_2001.log; [ -f $f ] || f=$f.gz
zcat -f $f | grep -a -E "Ast =|reco_sim: event" | head -4
