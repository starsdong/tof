#!/bin/bash
# usage: bash scripts/draw_qa.sh <inRoot> <plotDir>
IN=${1:-output/qa_test_10k.root}
DIR=${2:-plots/qa_test_10k}
cd /gpfs01/star/pwg/dongx/tof/TOFME || exit 1
bash scripts/sl7.sh "root -b -q -l 'macros/drawTofMEQA.C+(\"${IN}\",\"${DIR}\")'"
