#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
bash scripts/sl7.sh "cons" | grep -i -E " error|rror:|Install.*libStTofMEMaker.so"
mkdir -p logs/me_v9 output/me_v9
timeout 1500 bash scripts/run_me.sh 7000 filelists/me_v3/job_20091008.list output/me_v9_test.root logs/me_v9/test.log geom/btof_geom_v3.txt 5
grep -E "Finish|rror" logs/me_v9/test.log | tail -3
