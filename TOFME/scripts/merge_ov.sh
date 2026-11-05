#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
grep -a -h "StTofMEOverlayMaker::Finish" logs/mudst/ov_2*.log | cut -c1-250
bash scripts/sl7.sh "hadd -f output/mudst/ov.root output/mudst/ov_20060065_raw_1000011.root output/mudst/ov_20060065_raw_5500010.root output/mudst/ov_20073019_raw_4500013.root output/mudst/ov_20073019_raw_7000013.root" | tail -1
