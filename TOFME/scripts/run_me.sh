#!/bin/bash
# usage: bash scripts/run_me.sh <nEvents> <list> <outRoot> <log> [geomFile] [nMix]
# Runs macros/runTofME.C with root4star inside the SL7 container (SL24y).
N=${1:-10000}
LIST=${2:-filelists/test.list}
OUT=${3:-output/me_test.root}
LOG=${4:-logs/me_test.log}
GEOM=${5:-}
NMIX=${6:-5}
PDF=${7:-}
cd /gpfs01/star/pwg/dongx/tof/TOFME || exit 1
MACRO="macros/runTofME.C(${N},\"${LIST}\",\"${OUT}\",\"${GEOM}\",${NMIX},\"${PDF}\")"
echo "$MACRO" > "$LOG.cmd"
bash scripts/sl7.sh "root4star -b -q -l '${MACRO}' >&! ${LOG}"
echo "exit: $?"
