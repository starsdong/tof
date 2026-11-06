#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
ls output/sim/h50k/s*/tfs/simana.root > sim/h50k_ana.list; wc -l < sim/h50k_ana.list
bash scripts/sl7.sh "hadd -f -k output/sim/h50k_tfs_ana.root @sim/h50k_ana.list" | tail -1
ls -la output/sim/h50k_tfs_ana.root
du -sh output/sim/h50k 2>/dev/null | tail -1
ls output/sim/h1k/s*/trs/simana.root > sim/h1k_trs_ana.list; bash scripts/sl7.sh "hadd -f -k output/sim/h1k_trs_ana.root @sim/h1k_trs_ana.list" | tail -1
