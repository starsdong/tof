#!/bin/bash
# Stage the test picoDst to GPFS (faster iterative tests) and run the geometry pass + solve.
cd /gpfs01/star/pwg/dongx/tof/TOFME || exit 1
SRC=$(cat filelists/test.list)
DST=data/$(basename "$SRC")
if [ ! -s "$DST" ]; then xrdcp -s "$SRC" "$DST" || exit 1; fi
ls -la "$DST"
echo "$PWD/$DST" > filelists/test_gpfs.list
date +%T
bash scripts/run_me.sh 10000 filelists/test_gpfs.list output/me_geom_test.root logs/me_geom_test.log ""
date +%T
grep -E "runTofME|Finish|Error|error" logs/me_geom_test.log | grep -v GARFIELD | head
bash scripts/sl7.sh "root -b -q -l 'macros/makeGeom.C+(\"output/me_geom_test.root\",\"geom/btof_geom_test.txt\")'" 2>&1 | grep -v -E "GARFIELD|^$|Info in"
