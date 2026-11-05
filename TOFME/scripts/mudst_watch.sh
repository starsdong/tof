#!/bin/bash
# Background watcher on starsub01: wait for the Data Carousel restore of the test MuDst files, then run the probe.
# Start: nohup bash scripts/mudst_watch.sh > logs/mudst/watch.log 2>&1 &
# Checks every 30 min for up to 72 h. A file counts as complete when its size is unchanged over 3 min.
# Output: logs/mudst/watch.log, output/mudst_probe.root, marker logs/mudst/probe.done (contains OK or FAIL).
cd /gpfs01/star/pwg/dongx/tof/TOFME || exit 1
mkdir -p logs/mudst output
rm -f logs/mudst/probe.done
echo $$ > logs/mudst/watch.pid
for i in $(seq 1 144); do
  stable=""
  for f in mudst/*.MuDst.root; do
    [ -f "$f" ] || continue
    s1=$(stat -c %s "$f"); sleep 180; s2=$(stat -c %s "$f")
    if [ "$s1" == "$s2" ] && [ "$s2" -gt 1000000 ]; then stable="$f"; break; fi
  done
  echo "$(date '+%F %T') check $i: files $(ls mudst/*.MuDst.root 2>/dev/null | wc -l)/6, stable: ${stable:-none}"
  if [ -n "$stable" ]; then
    echo "$PWD/$stable" > filelists/mudst/test1.list
    echo "$(date '+%F %T') running probe on $stable"
    bash scripts/sl7.sh "root4star -b -q -l 'macros/probeMuBTof.C(5000,\"filelists/mudst/test1.list\",\"output/mudst_probe.root\")' >&! logs/mudst/probe.log"
    if grep -q "probeMuBTof: done" logs/mudst/probe.log; then st=OK; else st=FAIL; fi
    grep -E "probeMuBTof" logs/mudst/probe.log | tail -3
    echo "$(date '+%F %T') probe $st"
    echo "$st $stable $(date '+%F %T')" > logs/mudst/probe.done
    exit 0
  fi
  sleep 1800
done
echo "TIMEOUT no restored file after 72 h $(date '+%F %T')" > logs/mudst/probe.done
