#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
for m in 1 2 3; do grep -a -h "StTofMEOverlayMaker::Finish" logs/mudst/ov${m}_2*.log | cut -c1-230; bash scripts/sl7.sh "hadd -f output/mudst/ov${m}.root $(ls output/mudst/ov${m}_2*.root | tr "\n" " ")" | tail -1; done
