#!/bin/bash
# prepare the large MuDst production: per-file lists and tag files from filelists/mudst/restore_100.txt
cd /gpfs01/star/pwg/dongx/tof/TOFME
: > filelists/mudst/tags_big.txt; : > filelists/mudst/tags_big_ov.txt
for f in $(awk '{print $2}' filelists/mudst/restore_100.txt); do
  t=$(basename $f .MuDst.root); t=${t#st_physics_}
  echo "$f" > filelists/mudst/f_${t}.list
  echo "$t" >> filelists/mudst/tags_big.txt
  echo "$t, 1" >> filelists/mudst/tags_big_ov.txt
done
awk -F_ '{print $1}' filelists/mudst/tags_big.txt | sort -u > filelists/mudst/runs_big.txt
wc -l filelists/mudst/tags_big.txt filelists/mudst/runs_big.txt
