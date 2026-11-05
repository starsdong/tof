#!/bin/bash
# short interactive test of StTofMEMaker v2 (m2) on the 10k test file
cd /gpfs01/star/pwg/dongx/tof/TOFME || exit 1
N=${1:-3000}
bash scripts/run_me.sh $N filelists/test.list output/me_m2_test.root logs/me_m2_test.log geom/btof_geom_1M.txt 5
grep -E "Finish events|rror|StRefMultCorr bad|segmentation|Abort" logs/me_m2_test.log | grep -v -E "GARFIELD|Error Codes" | head -20
