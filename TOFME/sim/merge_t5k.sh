#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
for r in notot tot slow; do ls output/sim/t5k/s*/$r/simana.root > sim/t5k_$r.list; bash scripts/sl7.sh "hadd -f -k output/sim/t5k_${r}_ana.root @sim/t5k_$r.list" | tail -1; done
ls -la output/sim/t5k_*_ana.root
