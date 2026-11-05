#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
ls output/mudst/ | grep probe
bash scripts/sl7.sh "hadd -f output/mudst/probe4.root output/mudst/probe_20060065_raw_1000011.root output/mudst/probe_20060065_raw_5500010.root output/mudst/probe_20073019_raw_4500013.root output/mudst/probe_20073019_raw_7000013.root" | tail -1
