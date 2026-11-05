#!/bin/bash
# Build the TOFME QA note. Independent of the caller's cwd.
set -e
NOTE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$NOTE"
mkdir -p build tables logs figures
echo "note root: $NOTE"
if [ "$1" != "--tex-only" ]; then
  root -b -q -l "$NOTE/scripts/make_numbers.C(\"$NOTE/data/qa_test_10k.root\")" > logs/make_numbers.out 2>&1
  root -b -q -l "$NOTE/scripts/make_figures.C(\"$NOTE/data/qa_test_10k.root\")" > logs/make_figures.out 2>&1
  root -b -q -l "$NOTE/scripts/make_me.C(\"$NOTE/data/me_test.root\",\"$NOTE/data/btof_geom_test_qa.root\",\"test\")" > logs/make_me_test.out 2>&1
  root -b -q -l "$NOTE/scripts/make_me.C(\"$NOTE/data/me_v3.root\",\"$NOTE/data/btof_geom_v3_qa.root\",\"v3\")" > logs/make_me_v3.out 2>&1
  root -b -q -l "$NOTE/scripts/make_me_runs.C(\"$NOTE/data/me_v3_merged_inputs.txt\",\"$NOTE/data/me_v3_jobs\",\"$NOTE/data/btof_geom_v3.txt\")" > logs/make_me_runs.out 2>&1
  root -b -q -l "$NOTE/scripts/make_m2.C(\"$NOTE/data/me_v3.root\",\"v3\")" > logs/make_m2_v3.out 2>&1
  root -b -q -l "$NOTE/scripts/make_m2.C(\"$NOTE/data/me_v7.root\",\"v7\",\"\",true)" > logs/make_m2_v7.out 2>&1
  root -b -q -l "$NOTE/scripts/make_m2.C(\"$NOTE/data/me_v8.root\",\"v8\",\"\",true)" > logs/make_m2_v8.out 2>&1
  root -b -q -l "$NOTE/scripts/make_bkg.C(\"$NOTE/data/me_v4.root\")" > logs/make_bkg.out 2>&1
  root -b -q -l "$NOTE/scripts/make_cell.C(\"$NOTE/data/me_v5.root\")" > logs/make_cell.out 2>&1
  root -b -q -l "$NOTE/scripts/make_nb.C(\"$NOTE/data/me_v6.root\")" > logs/make_nb.out 2>&1
  root -b -q -l "$NOTE/scripts/make_late.C(\"$NOTE/data/me_v8.root\")" > logs/make_late.out 2>&1
  root -b -q -l "$NOTE/scripts/make_mu.C(\"$NOTE/data/me_mu.root\")" > logs/make_mu.out 2>&1
  root -b -q -l "$NOTE/scripts/make_ov.C(\"$NOTE/data/ov.root\")" > logs/make_ov.out 2>&1
  root -b -q -l "$NOTE/scripts/make_ov.C(\"$NOTE/data/ov1.root\",\"trk\",\"ovt\")" > logs/make_ovtrk.out 2>&1
  root -b -q -l "$NOTE/scripts/make_ov.C(\"$NOTE/data/ov2.root\",\"pool\",\"pool\")" > logs/make_ovpool.out 2>&1
  root -b -q -l "$NOTE/scripts/make_ov.C(\"$NOTE/data/ov3.root\",\"pooltrk\",\"poolt\")" > logs/make_ovpooltrk.out 2>&1
  root -b -q -l "$NOTE/scripts/make_probe.C(\"$NOTE/data/mudst_probe4.root\",\"$NOTE/data/closure_s0.root\")" > logs/make_probe.out 2>&1
  python3 "$NOTE/scripts/make_mudst_table.py" "$NOTE" > logs/make_mudst_table.out
  root -b -q -l "$NOTE/scripts/make_mu.C(\"$NOTE/data/mu_big.root\",\"big\")" > logs/make_mubig.out 2>&1
  root -b -q -l "$NOTE/scripts/make_qm2.C(\"$NOTE/data/mu_big.root\",1,\"mubig\",\"Emulation\")" > logs/make_qm2_mubig.out 2>&1
  for x in "mubig 0" "mubigneg -1" "mubigpos 1"; do set -- $x
    root -b -q -l "$NOTE/scripts/make_fit.C(\"$NOTE/data/mu_big.root\",1,\"$1\",\"Emul\",$2)" > logs/make_fit_$1.out 2>&1 &
  done; wait
  python3 "$NOTE/scripts/make_big_tables.py" "$NOTE" > logs/make_big_tables.out
  root -b -q -l "$NOTE/scripts/make_ov.C(\"$NOTE/data/ov_big.root\",\"big\",\"ovbig\")" > logs/make_ovbig.out 2>&1
  root -b -q -l "$NOTE/scripts/make_qm2.C(\"$NOTE/data/ov_big.root\",0,\"ovbig\",\"Overlay\")" > logs/make_qm2_ovbig.out 2>&1
  for x in "ovbig 0" "ovbigneg -1" "ovbigpos 1"; do set -- $x
    root -b -q -l "$NOTE/scripts/make_fit.C(\"$NOTE/data/ov_big.root\",0,\"$1\",\"Ovl\",$2)" > logs/make_fit_$1.out 2>&1 &
  done; wait
  (cd "$NOTE" && root -b -q -l scripts/make_lowp.C) > logs/make_lowp.out 2>&1
  python3 "$NOTE/scripts/make_big_m2_tables.py" "$NOTE" > logs/make_big_m2_tables.out
fi
(cd "$NOTE/.." && md5 -r StRoot/StTofMEQAMaker/StTofMEQAMaker.h StRoot/StTofMEQAMaker/StTofMEQAMaker.cxx macros/runTofMEQA.C macros/drawTofMEQA.C) > logs/md5_local.txt
python3 "$NOTE/scripts/parse_log.py" "$NOTE" > logs/parse_log.out
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex > logs/pdflatex_pass$i.out 2>&1 \
    || { echo "pdflatex pass $i FAILED"; tail -30 logs/pdflatex_pass$i.out; exit 1; }
done
echo "--- checks"
grep -c "Overfull" build/main.log | sed 's/^/overfull boxes: /' || true
grep -E "undefined|Missing|not found" build/main.log | grep -v "^Package: " | head || true
USED=$(grep -oh '\\figpanel{[0-9.]*}{[A-Za-z0-9_]*}' sections/*.tex | sed 's/.*}{//; s/}$//' | sort -u)
for b in $USED; do [ -f "figures/$b.pdf" ] || echo "MISSING figures/$b.pdf"; done
echo "figures used: $(echo $USED | wc -w)"
echo "unused figures:"
for f in figures/*.pdf; do b=$(basename "$f" .pdf); echo "$USED" | grep -qx "$b" || echo "  $b"; done
pdfinfo build/main.pdf 2>/dev/null | grep Pages || python3 -c "print('pages: n/a')"
echo "PDF: $NOTE/build/main.pdf"
