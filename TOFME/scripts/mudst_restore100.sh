#!/bin/bash
# Data Carousel request for ~100 MuDst files (P24iy, production_19GeV_2019): 7 files per run for the 14 runs of the
# picoDst v3 sample except 20088015 (trigger-data decoding problem). The raw files used in the picoDst sample
# (filelists/me_v3_files.txt) come first, then other files of the same run (largest first). Files already
# restored in mudst/ are skipped.  usage: bash scripts/mudst_restore100.sh [submit]
cd /gpfs01/star/pwg/dongx/tof/TOFME
OUT=filelists/mudst/restore_100.txt
: > $OUT
for r in $(cat filelists/me_v3_jobids.txt); do
  [ "$r" == "20088015" ] && continue
  get_file_list.pl -keys "path,filename,events" -cond "filetype=daq_reco_MuDst,production=P24iy,runnumber=$r,storage=HPSS" -limit 0 -delim " " 2>/dev/null | sort -k3 -n -r > /tmp/dongx_mu_$r.txt
  sel=0
  # files of the picoDst sample first
  for pico in $(grep "^$r " filelists/me_v3_files.txt | awk '{print $2}' | xargs -n1 basename); do
    mu=${pico/.picoDst.root/.MuDst.root}
    line=$(grep " $mu " /tmp/dongx_mu_$r.txt | head -1)
    [ -z "$line" ] && continue
    [ -f mudst/$mu ] && continue
    grep -q "/$mu " $OUT && continue
    [ $sel -ge 7 ] && break
    echo "$(echo $line | awk '{print $1"/"$2}') /gpfs01/star/pwg/dongx/tof/TOFME/mudst/$mu" >> $OUT; sel=$((sel+1))
  done
  while read p f e; do
    [ $sel -ge 7 ] && break
    [ -f mudst/$f ] && continue
    grep -q "/$f " $OUT && continue
    echo "$p/$f /gpfs01/star/pwg/dongx/tof/TOFME/mudst/$f" >> $OUT; sel=$((sel+1))
  done < /tmp/dongx_mu_$r.txt
  echo "run $r: $sel files"
  rm -f /tmp/dongx_mu_$r.txt
done
echo "total files: $(wc -l < $OUT)"
echo "events (catalog): $(awk '{print $2}' $OUT | xargs -n1 basename | while read f; do get_file_list.pl -keys events -cond "filetype=daq_reco_MuDst,production=P24iy,filename=$f,storage=HPSS" -limit 1 2>/dev/null; done | awk '{s+=$1} END {print s}')"
df -h /gpfs01/star/pwg/dongx | tail -1
if [ "$1" == "submit" ]; then hpss_user.pl -f $OUT; fi
