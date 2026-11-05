#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
bash scripts/sl7.sh "root4star -b -q -l 'macros/runTofMEMu.C(2000,\"filelists/mudst/run_20060065.list\",\"output/mudst/me_test.root\",\"geom/btof_geom_v3.txt\",5)' >&! logs/mudst/me_test.log"
grep -a -E "StTofMEMuMaker|runTofMEMu|rror|Make returned" logs/mudst/me_test.log | grep -v -E "GARFIELD|ShowMembers" | tail -8
