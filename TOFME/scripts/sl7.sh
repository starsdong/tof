#!/bin/bash
# Run a csh command inside the STAR SL7 container with library $STARLIB (default SL24y).
# usage: bash scripts/sl7.sh 'cons'      (run from /gpfs01/star/pwg/dongx/tof/TOFME)
STARLIB=${STARLIB:-SL24y}
SIF=/cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif
WORK=/gpfs01/star/pwg/dongx/tof/TOFME
cd "$WORK" || exit 1
exec apptainer exec -B /gpfs01,/gpfs,/star,/cvmfs,/direct,/sdcc "$SIF" \
  tcsh -c "starver $STARLIB >& /dev/null; cd $WORK; $1"
