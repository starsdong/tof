#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
bash scripts/sl7.sh "cons" | grep -i -E " error|Install.*libStTofSimAna"
ls output/sim/h1k/s*/tfs/hijing_*.MuDst.root > sim/h1k_tfs_mudst.list; wc -l < sim/h1k_tfs_mudst.list
bash scripts/sl7.sh "root4star -b -q -l 'macros/runTofSimAna.C(1000000,\"sim/h1k_tfs_mudst.list\",\"output/sim/h1k_tfs_ana.root\")' >&! logs/sim/h1k_tfs_ana.log"
grep -a "StTofSimAnaMaker::Finish" logs/sim/h1k_tfs_ana.log
