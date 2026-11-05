#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
N=${1:-300}; L=${2:-filelists/mudst/test1.list}; O=${3:-output/mudst/ov_test.root}; G=${4:-logs/mudst/ov_test.log}; M=${5:-0}
bash scripts/sl7.sh "root4star -b -q -l 'macros/runTofMEOverlay.C(${N},\"${L}\",\"${O}\",\"y2019\",0,${M})' >&! ${G}"
grep -a -E "StTofMEOverlayMaker v2|StTofMEOverlayMaker::Finish|QAInfo:Maker.*(btofMatch|tofMEov|btofCalib)|rror|Fatal|segmentation|break" ${G} | grep -v -E "GARFIELD|ShowMembers" | tail -8
