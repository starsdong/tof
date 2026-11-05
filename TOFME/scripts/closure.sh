#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
N=${1:-300}; S=${2:-0}
bash scripts/sl7.sh "root4star -b -q -l 'macros/runTofClosure.C(${N},\"filelists/mudst/test1.list\",\"output/mudst/closure_s${S}.root\",\"y2019\",${S})' >&! logs/mudst/closure_s${S}.log"
grep -a -E "StTofClosureMaker:|runTofClosure|rror|Fatal|tStart|VpdStart|startless|Start-less" logs/mudst/closure_s${S}.log | grep -v -E "GARFIELD|ShowMembers" | sort | uniq -c | sort -rn | head -25
