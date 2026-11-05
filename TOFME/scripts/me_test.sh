#!/bin/bash
# 10k test of the mixed-event setup: (1) geometry pass, (2) solve geometry, (3) projection + SE/ME pass.
# run detached:  nohup bash scripts/me_test.sh >& logs/me_test_driver.out < /dev/null &
cd /gpfs01/star/pwg/dongx/tof/TOFME || exit 1
LIST=${LIST:-filelists/test.list}
echo "start $(date +%T)"
if [ ! -s output/me_geom_test.root ]; then
  bash scripts/run_me.sh 10000 $LIST output/me_geom_test.root logs/me_geom_test.log ""
fi
echo "geom pass done $(date +%T)"
bash scripts/sl7.sh "root -b -q -l 'macros/makeGeom.C+(\"output/me_geom_test.root\",\"geom/btof_geom_test.txt\")' >&! logs/makeGeom_test.log"
grep -v -E "GARFIELD|^$" logs/makeGeom_test.log | tail -8
echo "geom solved $(date +%T)"
bash scripts/run_me.sh 10000 $LIST output/me_test.root logs/me_test.log geom/btof_geom_test.txt 5
echo "ME pass done $(date +%T)"
grep -E "runTofME|StTofMEMaker|rror" logs/me_test.log | grep -v GARFIELD | head -20
echo "END"
