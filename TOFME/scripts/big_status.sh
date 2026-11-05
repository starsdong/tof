#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
date
condor_q dongx -totals | grep query
np=0; no=0; nm=0
for t in $(cat filelists/mudst/tags_big.txt); do grep -a -q "probeMuBTof: done" logs/mudst/probe_${t}_p1.log 2>/dev/null && np=$((np+1)); done
for t in $(cat filelists/mudst/tags_all.txt); do
  grep -a -q "StTofMEOverlayMaker::Finish" logs/mudst/ov1_${t}.log 2>/dev/null && no=$((no+1))
  grep -a -q "StTofMEMuMaker::Finish" logs/mudst/mu_${t}.log 2>/dev/null && nm=$((nm+1))
done
echo "probe pass1 done $np/98 | method 2 done $no/102 | method 1 done $nm/102"
