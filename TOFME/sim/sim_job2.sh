#!/bin/bash
# sim_job2.sh seed nev (= sim_job.sh + KEEPFZD=0 deletes the .fzd; separate file so running sim_job.sh jobs are not edited under bash) [bmin bmax] : one HTCondor job of the HIJING Au+Au 19.6 GeV simulation
#  1. starsim (GEANT3, y2019a)            -> output/sim/h1k/s<seed>/gen/hijing_<seed>.fzd
#  2. bfc with the fast TPC simulator (tfs) -> .../tfs/hijing_<seed>.MuDst.root
#  3. bfc with the full TPC response (TpcRS) -> .../trs/hijing_<seed>.MuDst.root   (skipped if RECO=tfs)
#  4. StTofSimAnaMaker on each MuDst        -> .../{tfs,trs}/simana.root
SEED=$1; NEV=$2; BMIN=${3:-0}; BMAX=${4:-14.5}; RECO=${RECO:-both}; TOP=${TOP:-h1k}
WD=/gpfs01/star/pwg/dongx/tof/TOFME
BASE=$WD/output/sim/$TOP/s$SEED
COMMON="y2019a,fzin,AgML,useXgeom,iTpcIT,MakeEvent,ITTF,NoSsdIt,NoSvtIt,Idst,BAna,VFMinuit,tpcDB,btofsim,vpdSim,btofMatch,btofCalib,vpdCalib,IdTruth,McEvent,CMuDst,-hitfilt"
echo "sim_job: seed $SEED nev $NEV b $BMIN-$BMAX reco $RECO start $(date) on $(hostname)"
bash $WD/sim/run_sim.sh $NEV $SEED y2019a $BMIN $BMAX $BASE/gen gen
for r in tfs trs; do
  [ "$RECO" == "tfs" ] && [ $r == trs ] && continue
  mkdir -p $BASE/$r; ln -sf $BASE/gen/hijing_$SEED.fzd $BASE/$r/hijing_$SEED.fzd
  if [ $r == tfs ]; then CH="$COMMON,tfs"; else CH="$COMMON,TpcRS,TpxClu,TpcHitMover"; fi
  bash $WD/sim/run_sim.sh $NEV $SEED y2019a $BMIN $BMAX $BASE/$r reco "$CH"
  bash $WD/scripts/sl7.sh "root4star -b -q -l 'macros/runTofSimAna.C(1000000,\"$BASE/$r/hijing_$SEED.MuDst.root\",\"$BASE/$r/simana.root\")' >&! $BASE/$r/ana.log"
  grep -a "StTofSimAnaMaker::Finish" $BASE/$r/ana.log
  rm -f $BASE/$r/hijing_$SEED.event.root    # StEvent copy not needed
  gzip -f $BASE/$r/reco_$SEED.log
done
gzip -f $BASE/gen/sim_$SEED.log
# the .fzd (~3 MB/event) is only needed for re-reconstruction; delete it for large tfs-only productions
[ "${KEEPFZD:-1}" == "0" ] && rm -f $BASE/gen/hijing_$SEED.fzd $BASE/*/hijing_$SEED.fzd
echo "sim_job: done $(date)"
