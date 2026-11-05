#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
hpss_user.pl -f filelists/mudst/restore_100.txt 2>&1 | tail -3
head -2 filelists/mudst/restore_100.txt
