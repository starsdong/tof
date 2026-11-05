#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
condor_q dongx -totals | grep "query"
for f in logs/me_m3/me_*.log; do n=$(grep "events in chain" $f | awk "{print \$2}"); l=$(grep "  event " $f | tail -1); echo "$f total=$n last:$l"; done
echo "finished: $(grep -l "This is the end" logs/me_m3/me_*.log 2>/dev/null | wc -l)"
