#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
bash scripts/sl7.sh "root4star -b -q -l 'macros/probeMuBTof.C(5000,\"filelists/mudst/test1.list\",\"output/mudst_probe.root\")' >&! logs/mudst/probe.log"
grep -a -E "probeMuBTof" logs/mudst/probe.log | tail -3
