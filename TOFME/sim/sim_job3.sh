#!/bin/bash
# sim_job3.sh seed nev : HIJING Au+Au 19.6 GeV, GEANT3 y2019a, then the same GEANT events reconstructed three times with the fast
# TPC simulator and different BTOF responses of the patched StBTofSimMaker:
#   notot : fast BTOF simulator, ToT = 0 (as in h50k)
#   tot   : fast BTOF simulator + ToT from the cell dE (TOFME_TOT_DE0)
#   slow  : slow (avalanche) BTOF simulator with threshold + ToT from the avalanche charge (TOFME_SLOW, TOFME_TOT_Q0)
# Outputs output/sim/$TOP/s<seed>/{gen,notot,tot,slow}; StTofSimAnaMaker -> .../<variant>/simana.root. The .fzd is deleted at the end.
SEED=$1; NEV=$2; TOP=${TOP:-t5k}
WD=/gpfs01/star/pwg/dongx/tof/TOFME
BASE=$WD/output/sim/$TOP/s$SEED
CH="y2019a,fzin,AgML,useXgeom,iTpcIT,MakeEvent,ITTF,NoSsdIt,NoSvtIt,Idst,BAna,VFMinuit,tpcDB,btofsim,vpdSim,btofMatch,btofCalib,vpdCalib,IdTruth,McEvent,CMuDst,-hitfilt,tfs"
echo "sim_job3: seed $SEED nev $NEV start $(date) on $(hostname)"
bash $WD/sim/run_sim.sh $NEV $SEED y2019a 0 14.5 $BASE/gen gen
for r in notot tot slow; do
  mkdir -p $BASE/$r; ln -sf $BASE/gen/hijing_$SEED.fzd $BASE/$r/hijing_$SEED.fzd
  unset TOFME_TOT_DE0 TOFME_TOT_SIG TOFME_SLOW TOFME_TOT_Q0
  [ $r == tot ] && export TOFME_TOT_DE0=0.047 TOFME_TOT_SIG=0
  [ $r == slow ] && export TOFME_SLOW=1 TOFME_TOT_Q0=0.0095 TOFME_TOT_SIG=0
  bash $WD/sim/run_sim.sh $NEV $SEED y2019a 0 14.5 $BASE/$r reco "$CH"
  zcat -f $BASE/$r/reco_$SEED.log | grep -a "TOFME patch" | head -1
  bash $WD/scripts/sl7.sh "root4star -b -q -l 'macros/runTofSimAna.C(1000000,\"$BASE/$r/hijing_$SEED.MuDst.root\",\"$BASE/$r/simana.root\")' >&! $BASE/$r/ana.log"
  grep -a "StTofSimAnaMaker::Finish" $BASE/$r/ana.log
  rm -f $BASE/$r/hijing_$SEED.event.root $BASE/$r/hijing_$SEED.fzd
  gzip -f $BASE/$r/reco_$SEED.log
done
gzip -f $BASE/gen/sim_$SEED.log
rm -f $BASE/gen/hijing_$SEED.fzd
echo "sim_job: done $(date)"
