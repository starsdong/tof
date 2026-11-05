#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
for m in 1 2 3; do bash scripts/ov_test.sh 300 filelists/mudst/test1.list output/mudst/ov_test$m.root logs/mudst/ov_test$m.log $m > logs/mudst/ov_test${m}.sum 2>&1 & done
wait
for m in 1 2 3; do echo "== mode $m"; cat logs/mudst/ov_test${m}.sum; done
