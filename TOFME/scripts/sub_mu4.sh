#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
bash scripts/sl7.sh "root -b -q -l -e 'TFile*F=TFile::Open(\"output/mudst/mu_20057003_raw_7500004.root\"); for (const char* n : {\"hQM2SE\",\"hQM2MEprodB\",\"hM2SEneg\",\"hM2MEprodBneg\"}) printf(\"%s %s\\n\", n, F->Get(n) ? \"present\" : \"MISSING\");'" | grep -E "present|MISSING"
ls -la --time-style=+%T .sl73_x8664_gcc485/lib/libStTofMEMuMaker.so
condor_submit scripts/mubig4.sub | tail -1
bash scripts/big_status.sh
