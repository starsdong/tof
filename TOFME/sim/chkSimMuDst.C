// chkSimMuDst.C : quick content check of a simulation MuDst (tracks, BTOF hits/matches, idTruth, MC arrays)
void chkSimMuDst(const char* fn, int nev = 3)
{
  gROOT->Macro("loadMuDst.C");
  StMuDstMaker* mu = new StMuDstMaker(0, 0, "", fn, "", 1000);
  TChain* ch = mu->chain();
  printf("entries %lld\n", ch->GetEntries());
  for (int i = 0; i < nev && i < ch->GetEntries(); i++) {
    mu->Make();
    StMuDst* d = mu->muDst();
    StMuEvent* ev = d->event();
    int nG = d->numberOfGlobalTracks(), nP = d->numberOfPrimaryTracks(), nH = d->numberOfBTofHit();
    int nMatch = 0, nTruthTrk = 0, nHitTruth = 0, nGood = 0;
    for (int j = 0; j < nG; j++) { StMuTrack* t = d->globalTracks(j); if (t->idTruth() > 0) nTruthTrk++;
      if (t->index2BTofHit() >= 0) { nMatch++; const StMuBTofPidTraits& pt = t->btofPidTraits();
        if (nGood < 4) printf("   trk p %.2f beta %.3f tof %.2f L %.1f idTruth %d qa %d hit idTruth %d\n", t->p().mag(), pt.beta(), pt.timeOfFlight(), pt.pathLength(), t->idTruth(), t->qaTruth(), d->btofHit(t->index2BTofHit())->idTruth()); nGood++; } }
    for (int j = 0; j < nH; j++) if (d->btofHit(j)->idTruth() > 0) nHitTruth++;
    TClonesArray* mcT = d->mcArray(1); TClonesArray* mcV = d->mcArray(0);
    printf("ev %d run %d refMult %d vz %.1f | glob %d prim %d (idTruth>0 %d) | btofHits %d (idTruth>0 %d) matched globals %d | MC trk %d vtx %d | tStart %.2f\n",
           i, ev->runNumber(), ev->refMult(), ev->primaryVertexPosition().z(), nG, nP, nTruthTrk, nH, nHitTruth, nMatch,
           mcT ? mcT->GetEntriesFast() : -1, mcV ? mcV->GetEntriesFast() : -1, d->btofHeader() ? d->btofHeader()->tStart() : -999.);
    if (i == 0 && mcT && mcT->GetEntriesFast() > 0) { mcT->At(0)->Print(); if (mcV && mcV->GetEntriesFast() > 1) mcV->At(1)->Print(); }
  }
}
