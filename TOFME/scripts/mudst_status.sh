#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
echo "restored files:"; ls -la mudst/ | grep MuDst
echo "expected: $(wc -l < filelists/mudst/restore_test.txt)"
hpss_user.pl -w 2>&1 | tail -8
echo "== my requests"; hpss_user.pl -l 2>&1 | grep -i -E "20060065|20073019|20088015|status|pending|request" | head -15
