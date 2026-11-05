#!/bin/bash
# Build XRootD file lists from the catalog dump (local storage, latest production P24iy).
cd /gpfs01/star/pwg/dongx/tof/TOFME/filelists
XRD=root://xrdstar.rcf.bnl.gov:1095/
awk -v x=$XRD '{print x $2 "/" $3}' local_P24iy_all.txt > P24iy_19GeV_all_xrd.list
# test list: latest run, 1 file with >10k events
sort -k1,1nr local_P24iy_all.txt | awk -v x=$XRD '$4>10000 {print x $2 "/" $3}' | head -1 > test.list
wc -l P24iy_19GeV_all_xrd.list test.list
cat test.list
