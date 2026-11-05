#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
mkdir -p logs/mudst
nohup setsid bash scripts/mudst_watch.sh > logs/mudst/watch.log 2>&1 < /dev/null &
sleep 5
echo "pid $(cat logs/mudst/watch.pid)"; ps -o pid,etime,cmd -p $(cat logs/mudst/watch.pid)
cat logs/mudst/watch.log
