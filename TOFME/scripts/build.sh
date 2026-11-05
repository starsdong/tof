#!/bin/bash
cd /gpfs01/star/pwg/dongx/tof/TOFME
bash scripts/sl7.sh "cons" | grep -i -E " error|rror:|Install.*lib(StTofMEMaker|StTofMEMuMaker).so"
