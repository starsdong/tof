#!/bin/bash
# MuDst replicas per storage for run 20091008, and an existence check of the local (XRootD) copies
cd /gpfs01/star/pwg/dongx/tof/TOFME
for st in local NFS; do
  n=$(get_file_list.pl -keys "path,filename" -cond "filetype=daq_reco_MuDst,production=P24iy,runnumber=20091008,storage=$st" -limit 0 -delim " " 2>/dev/null | wc -l)
  echo "storage $st: $n entries"
done
get_file_list.pl -keys "path,filename,events,node" -cond "filetype=daq_reco_MuDst,production=P24iy,runnumber=20091008,storage=local" -limit 0 -delim " " 2>/dev/null | sort -k3 -n -r | head -4 > filelists/mudst/local_20091008_top.txt
cat filelists/mudst/local_20091008_top.txt
while read p f e n; do
  xrdfs root://xrdstar.rcf.bnl.gov:1095 stat "$p/$f" 2>&1 | grep -E "Size|error|No such" | head -1
done < filelists/mudst/local_20091008_top.txt
echo "== NFS spot check (first 3)"
get_file_list.pl -keys "path,filename" -cond "filetype=daq_reco_MuDst,production=P24iy,runnumber=20091008,storage=NFS" -limit 3 -delim " " 2>/dev/null | while read p f; do ls "$p/$f" 2>&1 | tail -1; done
echo "== HPSS entries for the run"
get_file_list.pl -keys "path,filename" -cond "filetype=daq_reco_MuDst,production=P24iy,runnumber=20091008,storage=HPSS" -limit 0 -delim " " 2>/dev/null | wc -l
