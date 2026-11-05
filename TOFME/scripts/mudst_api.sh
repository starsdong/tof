#!/bin/bash
# API details for the MuDst BTOF probe + one test MuDst file of run 20091008
cd /gpfs01/star/pwg/dongx/tof/TOFME
S=/star/nfs4/AFS/star/packages/SL24y/StRoot
echo "== StBTofHeader"; grep -n -E "tStart|vpdVz|fiberHeader|const;" $S/StEvent/StBTofHeader.h | head -15
echo "== StMuBTofPidTraits"; grep -n -E "timeOfFlight|beta|matchFlag|yLocal|zLocal|pathLength|position" $S/StMuDSTMaker/COMMON/StMuBTofPidTraits.h | head -12
echo "== StMuTrack btof"; grep -n -E "btofPidTraits|index2Global|id\(\)" $S/StMuDSTMaker/COMMON/StMuTrack.h | head -6
echo "== StMuDst btofHeader"; grep -n -E "btofHeader" $S/StMuDSTMaker/COMMON/StMuDst.h | head -3
echo "== StMuBTofHit associatedTrackId/index2Global doc"; sed -n 30,60p $S/StMuDSTMaker/COMMON/StMuBTofHit.h | grep -v "^ *$"
echo "== files run 20091008 (disk)"
get_file_list.pl -keys "path,filename,events,node,storage" -cond "filetype=daq_reco_MuDst,production=P24iy,runnumber=20091008,storage!=HPSS" -limit 0 -delim " " 2>/dev/null | sort -k3 -n -r | head -6
