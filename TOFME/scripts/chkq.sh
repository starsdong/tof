#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
grep -a "StTofMEMuMaker::Finish" logs/mudst/muq_test.log | cut -c1-200
bash scripts/sl7.sh "root -b -q -l macros/chkq.C" | grep entries
