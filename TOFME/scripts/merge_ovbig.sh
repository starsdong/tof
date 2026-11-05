#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
ls output/mudst/ov1_2*.root | wc -l
for t in $(cat filelists/mudst/tags_all.txt); do grep -a -c StTofMEOverlayMaker::Finish logs/mudst/ov1_$t.log; done | sort | uniq -c
bash scripts/sl7.sh "hadd -f -k output/mudst/ov_big.root $(ls output/mudst/ov1_2*.root | tr "\n" " ")" | tail -1
ls -la output/mudst/ov_big.root
