#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
ls output/mudst/mu_*.root | grep -v test | wc -l
bash scripts/sl7.sh "hadd -f -k output/mudst/mu_big.root $(ls output/mudst/mu_2*.root | tr "\n" " ")" | tail -1
ls -la output/mudst/mu_big.root
