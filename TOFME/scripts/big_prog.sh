#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
for t in $(cat filelists/mudst/tags_all.txt); do l=logs/mudst/ov1_$t.log; d=$(grep -a "  event " $l | tail -1 | awk "{print \$2}"); n=$(grep -a "runTofMEOverlay: processing" $l | awk "{print \$NF}"); echo "$t ${d:-0} ${n:-0}"; done > /tmp/dongx_ovprog.txt
awk "{d+=\$2; n+=\$3; f=\$2/\$3; if(min==\"\"||f<min)min=f; if(f>max)max=f} END {printf \"events done %d of %d (%.0f%%), per job min %.0f%% max %.0f%%\\n\", d, n, 100*d/n, 100*min, 100*max}" /tmp/dongx_ovprog.txt
condor_q dongx -af JobStartDate | sort -n | head -1
condor_q dongx -af JobStartDate | sort -n | tail -1
