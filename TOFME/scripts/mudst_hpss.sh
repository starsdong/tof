#!/bin/bash
# HPSS restore tools and MuDst availability on disk (local/XRootD) for any 19.6 GeV 2019 production
cd /gpfs01/star/pwg/dongx/tof/TOFME
echo "== tools"; which hpss_user.pl htar hsi 2>&1 | grep -v "no hpss"
hpss_user.pl -h 2>&1 | head -25
echo "== local MuDst per production (count of entries)"
for pr in P21ic P23id P24iy; do
  n=$(get_file_list.pl -keys "filename" -cond "filetype=daq_reco_MuDst,trgsetupname=production_19GeV_2019,production=$pr,storage=local" -limit 200 -delim " " 2>/dev/null | wc -l)
  echo "$pr local MuDst entries (first 200 max): $n"
done
echo "== HPSS MuDst for the 15 runs (P24iy): files and events"
for r in $(head -15 filelists/me_v3_jobids.txt); do
  get_file_list.pl -keys "filename,events" -cond "filetype=daq_reco_MuDst,production=P24iy,runnumber=$r,storage=HPSS" -limit 0 -delim " " 2>/dev/null | awk -v r=$r '{n++; e+=$2} END {print "run", r, "files", n, "events", e}'
done
echo "== quota/space"; df -h /gpfs01/star/pwg/dongx 2>&1 | tail -1
