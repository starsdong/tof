#!/bin/bash
# run_sim.sh nev seed tag bmin bmax outdir [stage] [chain]
# HIJING Au+Au 19.6 GeV with starsim (GEANT3, geometry tag) + full reconstruction with bfc, in the SL7 container.
# stage: all (default) | gen (starsim only) | reco (bfc on an existing hijing_<seed>.fzd in outdir)
# The local TOFME libraries (patched StBTofSimMaker) are used through a .sl73_x8664_gcc485 symlink in outdir.
# Outputs in outdir: hijing_<seed>.{fzd,gener.root,MuDst.root,geant.root,event.root}, sim_<seed>.log, reco_<seed>.log
NEV=$1; SEED=$2; TAG=${3:-y2019a}; BMIN=${4:-0}; BMAX=${5:-14.5}; OUT=$6; STAGE=${7:-all}
CHAIN=${8:-"y2019a,fzin,AgML,useXgeom,TpcRS,TpxClu,iTpcIT,MakeEvent,ITTF,NoSsdIt,NoSvtIt,Idst,BAna,VFMinuit,tpcDB,TpcHitMover,btofsim,vpdSim,btofMatch,btofCalib,vpdCalib,IdTruth,McEvent,geantout,evout,CMuDst,-hitfilt"}
WD=/gpfs01/star/pwg/dongx/tof/TOFME
mkdir -p $OUT; cd $OUT || exit 1
[ -e .sl73_x8664_gcc485 ] || ln -s $WD/.sl73_x8664_gcc485 .sl73_x8664_gcc485
F=hijing_$SEED
cp $WD/sim/starsim_hijing.C $WD/sim/reco_sim.C .
echo "run_sim: start $(date) host $(hostname) nev $NEV seed $SEED tag $TAG b $BMIN-$BMAX stage $STAGE"
if [ "$STAGE" != "reco" ]; then
  bash $WD/scripts/sl7.sh "cd $OUT && root4star -b -q 'starsim_hijing.C($NEV,$SEED,\"$F\",\"$TAG\",$BMIN,$BMAX,50.)' >&! sim_$SEED.log"
  ls -la $F.fzd
  echo "run_sim: starsim done $(date)"
fi
if [ "$STAGE" != "gen" ]; then
  bash $WD/scripts/sl7.sh "cd $OUT && root4star -b -q 'reco_sim.C($NEV,\"$F.fzd\",\"$CHAIN\")' >&! reco_$SEED.log"
  grep -a "reco_sim: finished" reco_$SEED.log
  grep -a "QAInfo:Chain" reco_$SEED.log | tail -1
  echo "run_sim: reco done $(date)"
fi
ls -la $OUT
