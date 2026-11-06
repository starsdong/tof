#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
bash scripts/sl7.sh "hadd -f output/sim/h1k_tfs_ana.root output/sim/h1k/s*/tfs/simana.root" | tail -1
ls -la output/sim/h1k_tfs_ana.root
