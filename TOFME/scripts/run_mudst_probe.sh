#!/bin/bash
# Run the MuDst BTOF probe on the restored files (inside the SL7 container).
# usage: bash scripts/run_mudst_probe.sh [nEvents]
cd /gpfs01/star/pwg/dongx/tof/TOFME || exit 1
N=${1:-5000}
ls -1 /gpfs01/star/pwg/dongx/tof/TOFME/mudst/*.MuDst.root > filelists/mudst/test.list 2>/dev/null
echo "restored files: $(wc -l < filelists/mudst/test.list)"
[ -s filelists/mudst/test.list ] || { echo "no restored MuDst yet"; exit 1; }
head -1 filelists/mudst/test.list > filelists/mudst/test1.list
MACRO="macros/probeMuBTof.C(${N},\"filelists/mudst/test1.list\",\"output/mudst_probe.root\")"
bash scripts/sl7.sh "root4star -b -q -l '${MACRO}' >&! logs/mudst/probe.log"
grep -E "probeMuBTof|rror" logs/mudst/probe.log | grep -v GARFIELD | tail -5
