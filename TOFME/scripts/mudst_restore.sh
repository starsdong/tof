#!/bin/bash
# Build and submit a Data Carousel request for the MuDst test: 2 files from each of 3 runs, the same raw files
# as used in the picoDst sample (filelists/me_v3_files.txt). Target: /gpfs01/star/pwg/dongx/tof/TOFME/mudst/
cd /gpfs01/star/pwg/dongx/tof/TOFME
mkdir -p mudst filelists/mudst
OUT=filelists/mudst/restore_test.txt
: > $OUT
for r in 20060065 20073019 20088015; do
  grep "^$r " filelists/me_v3_files.txt | head -2 | while read run url ev; do
    pico=$(basename $url)
    mu=${pico/.picoDst.root/.MuDst.root}
    hp=$(get_file_list.pl -keys "path,filename" -cond "filetype=daq_reco_MuDst,production=P24iy,filename=$mu,storage=HPSS" -limit 1 -delim "/" 2>/dev/null)
    [ -n "$hp" ] && echo "$hp /gpfs01/star/pwg/dongx/tof/TOFME/mudst/$mu" >> $OUT
  done
done
cat $OUT
echo "files: $(wc -l < $OUT)"
if [ "$1" == "submit" ]; then hpss_user.pl -f $OUT; fi
