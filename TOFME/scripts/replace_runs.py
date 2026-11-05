#!/usr/bin/env python3
"""Replace stalled runs of a selection with the nearest eligible run (same criteria as select_runs.py:
not a StRefMultCorr bad run, >= minEv events in <= 8 largest files) not yet used.
usage: replace_runs.py <tag> <minEv> <run1> [run2 ...]
appends filelists/<tag>/job_<newrun>.list, adds lines to filelists/<tag>_files.txt,
writes filelists/<tag>_jobids_repl.txt with the new run numbers (to queue)."""
import collections, re, sys

tag, minEv = sys.argv[1], int(sys.argv[2])
stalled = [int(x) for x in sys.argv[3:]]
maxFiles = 8
base = "/gpfs01/star/pwg/dongx/tof/TOFME/filelists"
xrd = "root://xrdstar.rcf.bnl.gov:1095/"
badh = open("/star/nfs4/AFS/star/packages/SL24y/StRoot/StRefMultCorr/BadRun.h").read()
blk = badh[badh.index("badrun_refmult_2019["):]
bad = set(int(x) for x in re.findall(r"\b(20\d{6})\b", blk[blk.index("{"):blk.index("};")]))
runs = collections.defaultdict(list)
for line in open(f"{base}/local_P24iy_all.txt"):
    run, path, fn, ev, node = line.split()
    runs[int(run)].append((int(ev), f"{path}/{fn}"))
sel = {}
for r, fl in runs.items():
    if r in bad or not (20056032 <= r <= 20093036):
        continue
    tot, take = 0, []
    for e, f in sorted(fl, reverse=True)[:maxFiles]:
        take.append((e, f)); tot += e
        if tot >= minEv:
            break
    if tot >= minEv:
        sel[r] = take
used = set(int(l) for l in open(f"{base}/{tag}_jobids.txt"))
good = sorted(sel)
new = []
with open(f"{base}/{tag}_files.txt", "a") as out:
    for s in stalled:
        cands = sorted((abs(r - s), r) for r in good if r not in used and r not in new)
        r = cands[0][1]
        new.append(r)
        with open(f"{base}/{tag}/job_{r}.list", "w") as jl:
            for e, f in sel[r]:
                out.write(f"{r} {xrd}{f} {e}\n"); jl.write(f"{xrd}{f}\n")
        print(f"stalled {s} -> {r}: {len(sel[r])} files, {sum(e for e, _ in sel[r])} events")
open(f"{base}/{tag}_jobids_repl.txt", "w").write("".join(f"{r}\n" for r in new))
