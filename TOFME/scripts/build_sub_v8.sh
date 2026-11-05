#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
bash scripts/sl7.sh "cons" | grep -i -E "rror|Install.*libStTofMEMaker.so"
mkdir -p output/me_v8 logs/me_v8
condor_submit scripts/me_v8.sub | tail -1
