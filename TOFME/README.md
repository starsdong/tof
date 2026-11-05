# Project: TOFME — BTOF mixed-event matching (TPC tracks <-> BTOF hits)

## Overview
Set up a mixed-event calculation for BTOF matching between TPC tracks and BTOF hits
(random/fake-match background estimate). Data: STAR picoDst (StRoot/StPicoEvent,
https://github.com/star-bnl/star-sw). Step 1: basic picoDst QA (tracks, TPC-TOF matched info).

## Environment (SDCC)
- Node: starsub01 (ssh alias, ControlMaster; AlmaLinux 9.8, login shell tcsh)
- Workdir: /gpfs01/star/pwg/dongx/tof/TOFME
- STAR libs are SL7 builds -> run inside container
  /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif via scripts/sl7.sh 'cmd' (starver SL24y)
- Build: bash scripts/sl7.sh cons  -> .sl73_x8664_gcc485/lib
- tcsh noclobber: use '>&!' to overwrite logs inside sl7.sh commands.

## Data
- Run19 Au+Au 19.6 GeV, trgsetupname=production_19GeV_2019; productions in catalog:
  P21ic/SL21c, P23id/SL23d, P24iy/SL24y (latest, used).
- P24iy picoDst: HPSS 100560 files; local 99273 files (1.42e9 evts); NFS 100958 entries
  but sampled latest-run NFS paths do not exist (stale catalog).
- Read via XRootD: root://xrdstar.rcf.bnl.gov:1095//home/starlib/... (DOUBLE slash, absolute path required)
- filelists/local_P24iy_all.txt (catalog dump), P24iy_19GeV_all_xrd.list, test.list (1 file, run 20093036, 17041 evts)
- MB trigger ids 640001/11/21/31/41/51 (test file: 9989/10000 events have 640051)

## Code
- StRoot/StTofMEQAMaker: QA maker (event, track, dE/dx, 1/beta, m2 vs p, yLocal/zLocal,
  tray/module occupancy, TOF match efficiency num/den)
- macros/runTofMEQA.C (root4star; reads only Event/Track/BTofHit/BTofPidTraits branches)
- macros/drawTofMEQA.C (+ user's style.C/draw.C) -> plots/<tag>/*.png,pdf
- scripts/run_qa.sh N list out log ; scripts/draw_qa.sh in plotdir ; scripts/make_lists.sh

## Status (2026-11-03)
- Step 1 QA done on 10k test events; analysis note analysis_note/build/main.pdf (draft v0.1, local Mac).
- Step 2 mixed-event setup done: StTofMEMaker; 10k test + 1M sample (54/57 jobs, 19 runs) merged in output/me_1M.root.
  Geometry geom/btof_geom_test.txt (and _1M). Note v0.3 (28 pages).
- Rebuild note: bash analysis_note/build.sh
See session logs.

## Todo
- [x] Mixed-event matching (cell level) + 1M batch run
- [ ] Low-pT projection (energy loss); pT dependence of P_ME; random-match correction for yields/m2
