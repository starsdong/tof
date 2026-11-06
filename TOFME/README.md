# Project: TOFME — random TPC-track <-> BTOF-hit matches from mixed events

## Overview
Estimate the random (combinatorial) part of the STAR BTOF matching between TPC tracks and BTOF hits with a
mixed-event (ME) technique, and subtract it from the same-event (SE) m^2 spectra.
Data: Run19 Au+Au 19.6 GeV (production_19GeV_2019), production P24iy (SL24y); picoDst first, then MuDst.
Analysis note: `analysis_note/` (draft v0.7, 78 pages; `TOFME_analysis_note_v0.7.pdf` on GitHub).
GitHub: https://github.com/starsdong/tof, folder TOFME/ (code, note sources, figures, PDF; no data/logs/tables).

## Status (2026-11-06): simulation study added (note draft v0.8, Sec. 7.10, 86 pages)
All planned studies done; results are in note v0.7 (Secs. 7.7-7.9). Open items are listed below and in the note,
Sec. "Open issues".

## Key results
### picoDst (Secs. 7.1-7.6, StTofMEMaker v1-v10)
- picoDst stores BTofHits only partly (65% of hits untimed), so a faithful ME is not possible there; the studies
  established the production matching rule and the need for MuDst.
- Production matching rule (SL24y): candidates are the projected cell +-1; cells claimed by several tracks are
  dropped; one candidate is taken; for several, largest ToT (< 40 ns) wins, then closest local y.

### MuDst, 4 files (Sec. 7.7)
- Closure: re-running StBTofMatchMaker + StBTofCalibMaker on MuDst (processMuDst) reproduces production
  (16,451/16,451 matches same cell; tof median -14 ps, DB timestamp not pinned).
- Method 2 (StTofMEOverlayMaker, mode 1 = overlay partner hits AND partner global tracks, re-run production makers)
  is the baseline. Hits-only overlay overestimates random matches 2.3x; a separate hit pool measures
  P(random | no own hit) and is only a cross-check.
- Method 1 (StTofMEMuMaker v3/v4, production emulation with partner-track claims, own per-cell calibration,
  RMS 0.43 ns) is the cross-check.

### MuDst, 102 files, 1,858,212 events, 643,544 used for m^2 (Secs. 7.8, 7.9)
| | method 1 | method 2 |
|---|---|---|
| random matches 0-5% | 2.76% | 2.64% |
| random matches 30-40% | 0.84% | 0.80% |
| random matches 60-80% | 0.16% | 0.14% |
| random matches 0-80% | 1.65% | 1.58% |
| SE/ME between p and d (1.6<m^2<2.8), 0-5% | 1.32 | 1.41 |
- Under the peaks (method 2, 0-5%, 0.8-1 GeV/c): pi 1.7%, K 2.0%, p 0.9%; at 0.2-0.4 GeV/c K 11%, p 23%.
- SE-ME fits (Poisson likelihood of SE = signal + const + fixed normalized ME; G, double G, Student-t, split
  Student-t): the split Student-t is best in 8-9 of 9 momentum bins (proton high-m^2 tail). Mean chi2/ndf 0-5%:
  method 2 T2 12 (G 211). No model is statistically acceptable below 2 GeV/c (resolution mixing in each p bin).
- Yields of the two methods agree within 1.8% (pi, K, p, 0.4-2 GeV/c); raw ratios 0-5% at 0.4-2.5 GeV/c:
  pi-/pi+ 1.00-1.07, K-/K+ 0.60-0.62, pbar/p 0.092-0.104 (no efficiency correction).
- Low-momentum SE excess between p and d (Sec. 7.9.3): SE-ME per 10^4 tracks is flat in centrality
  (2.1-2.3 at 0.2-0.4 GeV/c, method 2) while the ME rises 28x -> not random matches; flat in t - t_pi up to 27 ns;
  charge-symmetric at 0.3-0.6 GeV/c (neg/pos 0.75-0.92), positive-only at 0.2-0.25 GeV/c (proton energy-loss tail).
  Interpretation (unproven): track-correlated late hits (secondaries/back-splash or delayed neighbour signals).

### Full simulation, HIJING+GEANT (Sec. 7.10, 2026-11-05/06)
- 49,857 HIJING Au+Au 19.6 GeV events (b 0-14.5 fm), GEANT3 y2019a, fast TPC simulator (tfs); 1k events also with TpcRS
  (category fractions agree within 1.6 sigma). Truth classes of the matched hit: own / decay / secondary / relative / unrelated.
- 0.2-0.4 GeV/c: decay muons 2.4% (stay near the pion peak), secondaries 2.7%, 100% of them starting in the BTOF trays
  (R 200-240 cm). In 1.6<m^2<2.8 the correlated rate is 2.7 per 1e4 QA tracks, flat in centrality; data/sim 0.71-1.00.
- Not reproduced: 0.4-0.6 GeV/c (data/sim 1.6-2.4), momentum-independent late plateau, charge pattern.
- Data mass spectrum vs ME + simulated correlated matches (Sec. 7.10.3): between p and d at 0.2-0.4 GeV/c SE/(ME+sim) =
  0.86-1.00 (SE/ME 1.75-10.9), but the simulation over-predicts the K-p valley and m2 > 3 by 2-3x -> proposed: evidence, not a
  correction (decision pending).
- ToT test (Sec. 7.10.4; 5k events, same GEANT events with no ToT / fast + ToT(dE) / slow + ToT(charge); patched
  StBTofSimMaker env vars TOFME_TOT_DE0, TOFME_SLOW, TOFME_TOT_Q0): ToT raises matching efficiency 2-3 points, doubles to
  triples unrelated matches and the 0.4-0.6 GeV/c correlated rate; it does not remove the K-p valley over-prediction.
  Code: sim/sim_job3.sh, sim/sim3.sub, macros/totData.C (data ToT), scripts/make_simtmpl.C, make_simvar.C.
- Code: sim/ (starsim_hijing.C, reco_sim.C, run_sim.sh, sim_job*.sh, sim*.sub, sim_status.sh), StRoot/StBTofSimMaker (patched:
  truth id for every hit), StRoot/StTofSimAnaMaker, macros/runTofSimAna.C, analysis_note/scripts/make_sim.C.
  Outputs on SDCC output/sim/{h1k,h50k}; merged analysis files analysis_note/data/sim_{h50k_tfs,h1k_tfs,h1k_trs}.root.

## Code
- StRoot/StTofMEQAMaker (QA), StTofMEMaker (picoDst ME, v10), StTofClosureMaker (MuDst closure),
  StTofMEMuMaker (method 1, v4), StTofMEOverlayMaker (method 2, v3; setMode 0 hits / 1 hits+tracks / 2 pool /
  3 pool+tracks). Both MuDst makers store hM2SE/ME, hQM2SE/ME (q*m^2) and hM2SEneg/MEneg (negative tracks).
- macros/: runTofMEQA.C, runTofMEMu.C, runTofMEOverlay.C (N, list, out, geom, 0, mode), runTofClosure.C,
  probeMuBTof.C (per-cell time offsets), chkq.C.
- scripts/ (SDCC, HTCondor): big_prep.sh (per-file lists), probe1.sub + off_big.sh (calibration pass -> off_<run>.root),
  mubig.sub (method 1), ovbig.sub (method 2), big_status.sh / big_prog.sh, merge_mubig.sh / merge_ovbig.sh,
  mudst_restore100.sh (HPSS restore of 98 MuDst files, list filelists/mudst/restore_100.txt).
- analysis_note/scripts/: make_mu.C (method 1, tag), make_ov.C (method 2, tag), make_qm2.C, make_fit.C
  (fn, method, tag, macro prefix, charge), make_lowp.C, make_big_tables.py, make_big_m2_tables.py, plus the
  earlier picoDst scripts. `bash analysis_note/build.sh` regenerates every table/figure and the PDF
  (`--tex-only` for LaTeX only). All numbers in the note are script-generated.

## Environment (SDCC)
- Node starsub01 (login shell tcsh: put loops/redirects in bash scripts piped via `cat >!`); workdir
  /gpfs01/star/pwg/dongx/tof/TOFME. STAR libs are SL7: run in the container via `bash scripts/sl7.sh 'cmd'`
  (starver SL24y); build `bash scripts/sl7.sh cons`.
- Local Mac copy ~/Work/Hermes/TOFME is the source of truth; push code with
  `COPYFILE_DISABLE=1 tar czf - --no-xattrs StRoot macros scripts | ssh starsub01 'cd /gpfs01/star/pwg/dongx/tof/TOFME && tar xzf -'`.
- Outputs: SDCC output/mudst/{mu_big,ov_big}.root; local analysis_note/data/{mu_big,ov_big}.root.
- picoDst/MuDst via XRootD root://xrdstar.rcf.bnl.gov:1095//home/starlib/... (double slash). Node qgp017 stalls; run
  20088015 MuDst has trigger-decoding floods (excluded).

## Open issues / possible next steps
- [ ] Identify the track-correlated late-hit source (cell of the late hit vs projection, ToT, |eta|); decide on a
      per-track template correction (e.g. from peripheral events).
- [ ] Line shapes: finer momentum bins or a resolution-scaled variable (1/beta residuals).
- [ ] Method 2 loses 3.6% of own matches in the overlay; quantify its effect on the random fraction.
- [ ] Pin the DB timestamp to the production snapshot (closure -14 ps).

Session history: session_log_2026-11-03.md (local only).
