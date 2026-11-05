#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
bash scripts/sl7.sh "cons" | grep -i -E " error|rror:|Install.*lib(StTofMEOverlayMaker|StTofMEMuMaker).so"
bash scripts/ov_test.sh 200 filelists/mudst/test1.list output/mudst/ovq_test.root logs/mudst/ovq_test.log 1 | grep -E "Finish|rror|break" &
bash scripts/sl7.sh "root4star -b -q -l 'macros/runTofMEMu.C(400,\"filelists/mudst/test1.list\",\"output/mudst/muq_test.root\",\"geom/btof_geom_v3.txt\",5)' >&! logs/mudst/muq_test.log"; grep -a -E "StTofMEMuMaker::Finish|rror|break" logs/mudst/muq_test.log | head -3
wait
bash scripts/sl7.sh "root -b -q -l -e 'for (const char* f : {\"output/mudst/ovq_test.root\",\"output/mudst/muq_test.root\"}) { TFile* F = TFile::Open(f); for (const char* n : {\"hQM2SE\",\"hQM2ME\",\"hQM2MEprodB\"}) { TH1* h = (TH1*)F->Get(n); if (h) printf(\"%s %s entries %.0f mean %.3f\\n\", f, n, h->GetEntries(), ((TH3*)h)->GetMean(3)); } }'" | grep entries
