#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
bash scripts/big_prep.sh
condor_submit scripts/ovbig.sub | tail -1
condor_submit scripts/probe1.sub | tail -1
