#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
date
L=filelists/mudst/restore_100.txt
n=0; s=0; for f in $(awk "{print \$2}" $L); do if [ -f $f ]; then n=$((n+1)); s=$((s+$(stat -c %s $f))); fi; done
echo "restored: $n / $(wc -l < $L) files, $((s/1024/1024/1024)) GB"
ls -lt --time-style=+%m-%d_%H:%M mudst/ | head -4
hpss_user.pl -w 2>&1 | tail -6
