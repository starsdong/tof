#!/bin/bash
# MuDst availability for production_19GeV_2019 P24iy, and the StMuBTofHit API in SL24y
cd /gpfs01/star/pwg/dongx/tof/TOFME
mkdir -p logs/mudst filelists/mudst
for st in local NFS HPSS; do
  echo "== storage $st"
  get_file_list.pl -keys "count(filename),sum(events)" -cond "filetype=daq_reco_MuDst,trgsetupname=production_19GeV_2019,production=P24iy,storage=$st" -limit 0 2>&1 | tail -2
done
echo "== runs of the v3 selection with MuDst on local/NFS"
for r in $(head -15 filelists/me_v3_jobids.txt); do
  n=$(get_file_list.pl -keys "path,filename,events,node,storage" -cond "filetype=daq_reco_MuDst,production=P24iy,runnumber=$r,storage!=HPSS" -limit 0 -delim " " 2>/dev/null | wc -l)
  echo "run $r disk MuDst files: $n"
done
echo "== StMuBTofHit.h"
S=/star/nfs4/AFS/star/packages/SL24y/StRoot
grep -n -E "^\s+(int|unsigned|float|double|Int_t|UInt_t|Float_t|Double_t|UChar_t|UShort_t|Short_t|bool|Bool_t)[^;]*\(\) *const" $S/StMuDSTMaker/COMMON/StMuBTofHit.h | head -30
grep -n -E "^\s+(Int_t|UInt_t|Float_t|Double_t|UChar_t|UShort_t|Short_t|Char_t) +m" $S/StMuDSTMaker/COMMON/StMuBTofHit.h
echo "== StMuDst btof accessors"
grep -n -i "btofhit\|numberOfBTofHit" $S/StMuDSTMaker/COMMON/StMuDst.h | head
