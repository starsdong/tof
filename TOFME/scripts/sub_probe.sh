#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
mkdir -p output/mudst
ls mudst/*.MuDst.root | sed -e "s|mudst/st_physics_||" -e "s|.MuDst.root||" > filelists/mudst/tags.txt
cat filelists/mudst/tags.txt
condor_submit scripts/probe.sub | tail -1
