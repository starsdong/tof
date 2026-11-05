#!/bin/bash
# Print software provenance for the TOFME QA (run on starsub01).
cd /gpfs01/star/pwg/dongx/tof/TOFME || exit 1
echo "host: $(hostname)  os: $(cat /etc/redhat-release)"
bash scripts/sl7.sh 'echo "container os: `cat /etc/redhat-release`"; echo "STAR=$STAR"; echo "STAR_LEVEL=$STAR_LEVEL"; echo "STAR_HOST_SYS(lib dir)=`ls -d .sl73_*`"; which root4star; echo "ROOTSYS=$ROOTSYS"; root-config --version; gcc --version | head -1; cat $STAR/.Version'
echo "container: /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif -> $(readlink -f /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif)"
md5sum StRoot/StTofMEQAMaker/* macros/runTofMEQA.C macros/drawTofMEQA.C
