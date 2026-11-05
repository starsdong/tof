#!/bin/bash
# usage: bash scripts/run_qa.sh <nEvents> <list> <outRoot> <log>
# Runs macros/runTofMEQA.C with root4star inside the SL7 container (SL24y).
N=${1:-10000}
LIST=${2:-filelists/test.list}
OUT=${3:-output/qa_test_10k.root}
LOG=${4:-logs/qa_test_10k.log}
cd /gpfs01/star/pwg/dongx/tof/TOFME || exit 1
MACRO="macros/runTofMEQA.C(${N},\"${LIST}\",\"${OUT}\")"
echo "$MACRO" > "$LOG.cmd"
bash scripts/sl7.sh "root4star -b -q -l '${MACRO}' >&! ${LOG}"
echo "exit: $?"
