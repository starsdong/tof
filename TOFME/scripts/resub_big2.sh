#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
condor_rm 817661
bash scripts/sl7.sh "cons" | grep -i -E " error|rror:|Install.*lib(StTofMEOverlayMaker|StTofMEMuMaker).so"
bash scripts/ov_test.sh 100 filelists/mudst/test1.list output/mudst/ovq_test.root logs/mudst/ovq_test.log 1 | grep -E "Finish|rror:|break"
bash scripts/sl7.sh "root4star -b -q -l 'macros/runTofMEMu.C(200,\"filelists/mudst/test1.list\",\"output/mudst/muq_test.root\",\"geom/btof_geom_v3.txt\",5)' >&! logs/mudst/muq_test.log"; grep -a "StTofMEMuMaker::Finish" logs/mudst/muq_test.log | cut -c1-120
for t in $(cat filelists/mudst/tags_all.txt); do rm -f output/mudst/ov1_${t}.root logs/mudst/ov1_${t}.log; done
condor_submit scripts/ovbig.sub | tail -1
