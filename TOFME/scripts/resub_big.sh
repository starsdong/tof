#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
condor_rm 817659
sleep 5
cat filelists/mudst/tags_ov.txt filelists/mudst/tags_big.txt > filelists/mudst/tags_all.txt
awk '{print $1", 1"}' filelists/mudst/tags_all.txt > filelists/mudst/tags_big_ov.txt
sed -i 's|queue tag from filelists/mudst/tags_big.txt|queue tag from filelists/mudst/tags_all.txt|' scripts/mubig.sub
for t in $(cat filelists/mudst/tags_all.txt); do rm -f output/mudst/ov1_${t}.root logs/mudst/ov1_${t}.log; done
wc -l filelists/mudst/tags_all.txt filelists/mudst/tags_big_ov.txt; grep queue scripts/mubig.sub scripts/ovbig.sub
condor_submit scripts/ovbig.sub | tail -1
